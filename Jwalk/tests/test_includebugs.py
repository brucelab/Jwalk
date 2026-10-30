#===============================================================================
#     This file is part of Jwalk.
#
#     Tests for -includebugs, which makes the optimized code reproduce the
#     output of the original Jwalk: the original path extraction (a SASD is
#     overestimated when two residues share a surface voxel) and the original
#     surface spheres chosen by the first letter of the atom name.
#
#     Run from the Jwalk directory with the launcher's virtual environment:
#
#         src/venv/bin/python -m unittest discover -s tests -v
#===============================================================================

import os
import shutil
import subprocess
import sys
import unittest

JWALK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(JWALK_DIR, "src")
EXAMPLES = os.path.join(JWALK_DIR, "Examples")
TEST_PDB = os.path.join(EXAMPLES, "test.pdb")
JWALK_SCRIPT = os.path.join(SRC_DIR, "Jwalk", "jwalk")

sys.path.insert(0, SRC_DIR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Jwalk import PDBTools, GridTools, SASDTools
from test_ligand import run_jwalk, read_crosslink_list


def read_paths(directory, pdb):
    """{pair: number of atoms in the path} from Jwalk_results/<pdb>_crosslinks.pdb"""
    name = os.path.splitext(os.path.basename(pdb))[0]
    paths = {}
    with open(os.path.join(directory, "Jwalk_results", name + "_crosslinks.pdb")) as f:
        for line in f:
            if line.startswith("MODEL"):
                pair = line.split()[2]
                paths[pair] = 0
            elif line.startswith("ATOM"):
                paths[pair] += 1
    return paths


class OriginalPaths(unittest.TestCase):
    """The one lysine pair of test.pdb whose SASD the original overestimated (performance report)."""

    @classmethod
    def setUpClass(cls):
        cls.fixed_dir, _ = run_jwalk(TEST_PDB)
        cls.bugs_dir, cls.bugs_out = run_jwalk(TEST_PDB, "-includebugs")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.fixed_dir)
        shutil.rmtree(cls.bugs_dir)

    def test_same_pairs(self):
        fixed = read_crosslink_list(self.fixed_dir, TEST_PDB)
        bugs = read_crosslink_list(self.bugs_dir, TEST_PDB)
        self.assertEqual(len(fixed), 66)
        self.assertEqual(set(fixed), set(bugs))
        self.assertIn("66 SASDs calculated", self.bugs_out)

    def test_only_125_135_differs_and_is_longer(self):
        fixed = read_crosslink_list(self.fixed_dir, TEST_PDB)
        bugs = read_crosslink_list(self.bugs_dir, TEST_PDB)
        changed = [k for k in fixed if abs(fixed[k] - bugs[k]) > 1e-6]
        self.assertEqual(changed, [("LYS-125-A-CA", "LYS-135-A-CA")])
        self.assertAlmostEqual(fixed[changed[0]], 15.948, places=3)
        self.assertAlmostEqual(bugs[changed[0]], 16.226, places=3)

    def test_duplicated_c_alpha_in_paths(self):
        fixed = read_paths(self.fixed_dir, TEST_PDB)
        bugs = read_paths(self.bugs_dir, TEST_PDB)
        self.assertEqual(set(fixed), set(bugs))
        longer = [p for p in fixed if bugs[p] > fixed[p]]
        self.assertTrue(longer)
        # a path that grew did so by whole C-alpha pairs
        self.assertTrue(all((bugs[p] - fixed[p]) % 2 == 0 for p in longer))
        self.assertFalse([p for p in fixed if bugs[p] < fixed[p]])

    def test_in_process_matches_command_line(self):
        structure = PDBTools.read_PDB_file(TEST_PDB)
        grid = GridTools.makeGrid(structure, 1)
        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "LYS", "LYS")
        dens_map = GridTools.generate_solvent_accessible_surface(grid, structure, aa1_CA, aa2_CA)
        aa1_voxels, _ = GridTools.find_surface_voxels(aa1_CA, dens_map, False, [])
        aa2_voxels, _ = GridTools.find_surface_voxels(aa2_CA, dens_map, False, [])
        sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, [], 60, 1, 1, [],
                                       original=True)
        sasds = GridTools.remove_duplicates(sasds)
        by_pair = {(k[0][0], k[1][0]): k[2] for k in sasds}
        self.assertAlmostEqual(by_pair[(125, 135)], 16.226, places=3)
        # two cpus give the same as one
        sasds2 = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, [], 60, 1, 2, [],
                                        original=True)
        self.assertEqual(sasds, GridTools.remove_duplicates(sasds2))


class LegacySpheres(unittest.TestCase):
    """With legacy_spheres only C, N, O and S atoms (by name) are placed in the surface."""

    @classmethod
    def setUpClass(cls):
        cls.base = PDBTools.read_PDB_file(TEST_PDB)
        extreme = cls.base.get_extreme_values()
        cls.far = (extreme[1] + 8, (extreme[2] + extreme[3]) / 2, (extreme[4] + extreme[5]) / 2)

    def voxels(self, structure, legacy):
        grid = GridTools.makeGrid(structure, 1)
        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "LYS", "LYS")
        dens_map = GridTools.generate_solvent_accessible_surface(grid, structure, aa1_CA, aa2_CA,
                                                                 legacy_spheres=legacy)
        return int((dens_map.fullMap > 0).sum())

    def added(self, elem, name, legacy):
        atom = self.base.atomList[0].copy()
        atom.elem, atom.atom_name, atom.res, atom.res_no, atom.chain = elem, name, "LIG", 999, "Z"
        atom.record_name = "HETATM"
        atom.x, atom.y, atom.z = self.far
        structure = PDBTools.BioPy_Structure(list(self.base.atomList) + [atom])
        return self.voxels(structure, legacy) - self.voxels(self.base, legacy)

    def test_protein_surface_unchanged(self):
        self.assertEqual(self.voxels(self.base, True), self.voxels(self.base, False))

    def test_protein_atoms_same_spheres(self):
        for elem, name in (("C", "C1"), ("N", "N1"), ("O", "O1"), ("S", "SD")):
            self.assertEqual(self.added(elem, name, True), self.added(elem, name, False), name)

    def test_other_atoms_get_no_sphere(self):
        for elem, name in (("P", "P"), ("FE", "FE1"), ("MO", "MO1"), ("ZN", "ZN"), ("MG", "MG")):
            self.assertEqual(self.added(elem, name, True), 0, name)
            self.assertGreater(self.added(elem, name, False), 0, name)

    def test_selenium_gets_sulfur_sphere_by_name(self):
        self.assertEqual(self.added("SE", "SE", True), 33)   # 2-voxel sphere, as S
        self.assertEqual(self.added("SE", "SE", False), 123)  # 3-voxel sphere, as Se


class CommandLine(unittest.TestCase):

    def test_help_lists_includebugs(self):
        result = subprocess.run([sys.executable, JWALK_SCRIPT, "-h"], env=dict(os.environ, PYTHONPATH=SRC_DIR),
                                capture_output=True, text=True)
        self.assertIn("-includebugs", result.stdout)


if __name__ == "__main__":
    unittest.main()
