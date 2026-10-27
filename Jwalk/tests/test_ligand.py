#===============================================================================
#     This file is part of Jwalk.
#
#     Tests for ligand mode (-ligand): reading hetero atoms, the element based
#     surface spheres, the freesasa .rsa parser, and a ligand placed on a
#     crosslink path in Examples/test_ligand.pdb.
#
#     Run from the Jwalk directory with the launcher's virtual environment:
#
#         src/venv/bin/python -m unittest discover -s tests -v
#
#     The -surface tests are skipped when freesasa is not on the PATH.
#===============================================================================

import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

JWALK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(JWALK_DIR, "src")
EXAMPLES = os.path.join(JWALK_DIR, "Examples")
TEST_PDB = os.path.join(EXAMPLES, "test.pdb")
LIGAND_PDB = os.path.join(EXAMPLES, "test_ligand.pdb")
JWALK_SCRIPT = os.path.join(SRC_DIR, "Jwalk", "jwalk")

sys.path.insert(0, SRC_DIR)

from Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools


def lysine_sasds(pdb, ligand, vox=1):
    """
    Runs the Jwalk pipeline in process on every lysine pair of pdb.

    Returns ({(res_no_1, res_no_2): sasd}, grid origin).
    """
    structure = PDBTools.read_PDB_file(pdb, hetatm=ligand)
    grid = GridTools.makeGrid(structure, vox)
    aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "LYS", "LYS")
    dens_map = GridTools.generate_solvent_accessible_surface(grid, structure, aa1_CA, aa2_CA)
    aa1_voxels, _ = GridTools.find_surface_voxels(aa1_CA, dens_map, False, [])
    aa2_voxels, _ = GridTools.find_surface_voxels(aa2_CA, dens_map, False, [])
    sasds = SASDTools.parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, [],
                                   60, vox, 1, [])
    sasds = GridTools.remove_duplicates(sasds)
    return {(k[0][0], k[1][0]): k[2] for k in sasds}, list(grid.origin)


def surface_voxels(structure):
    """Number of voxels inside the solvent accessible surface of structure, at 1 A."""
    grid = GridTools.makeGrid(structure, 1)
    aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "LYS", "LYS")
    dens_map = GridTools.generate_solvent_accessible_surface(grid, structure, aa1_CA, aa2_CA)
    return int((dens_map.fullMap > 0).sum())


def run_jwalk(pdb, *options):
    """
    Runs the jwalk command line on a copy of pdb in a fresh directory.

    Returns (directory, stdout). The directory holds Jwalk_results/.
    """
    directory = tempfile.mkdtemp(prefix="jwalk_test_")
    shutil.copy(pdb, directory)
    env = dict(os.environ, PYTHONPATH=SRC_DIR)
    result = subprocess.run([sys.executable, JWALK_SCRIPT, "-i", os.path.basename(pdb), "-ncpus", "1"]
                            + list(options),
                            cwd=directory, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError("jwalk failed:\n" + result.stdout + result.stderr)
    return directory, result.stdout


def read_crosslink_list(directory, pdb):
    """{(atom1, atom2): sasd} from Jwalk_results/<pdb>_crosslink_list.txt"""
    name = os.path.splitext(os.path.basename(pdb))[0]
    out = {}
    with open(os.path.join(directory, "Jwalk_results", name + "_crosslink_list.txt")) as f:
        next(f)
        for line in f:
            col = line.split()
            out[(col[2], col[3])] = float(col[4])
    return out


class ReadPDB(unittest.TestCase):

    def test_hetatm_off_by_default(self):
        structure = PDBTools.read_PDB_file(TEST_PDB)
        self.assertEqual(len(structure), 1001)
        self.assertTrue(all(a.record_name == "ATOM" for a in structure.atomList))

    def test_hetatm_reads_ligands_but_not_water(self):
        structure = PDBTools.read_PDB_file(TEST_PDB, hetatm=True)
        hetero = [a for a in structure.atomList if a.record_name == "HETATM"]
        self.assertEqual(len(structure), 1019)
        self.assertEqual(sorted(set(a.res for a in hetero)), ["BME", "SE4"])
        self.assertNotIn("HOH", set(a.res for a in structure.atomList))

    def test_element_is_read(self):
        structure = PDBTools.read_PDB_file(TEST_PDB, hetatm=True)
        elements = set(a.elem for a in structure.atomList)
        self.assertEqual(elements, {"C", "N", "O", "S", "SE"})


class ElementSpheres(unittest.TestCase):
    """Each element gets the expected sphere when placed in empty space."""

    SPHERE_VOXELS = {3: 123, 2: 33, 0: 0}

    @classmethod
    def setUpClass(cls):
        cls.base = PDBTools.read_PDB_file(TEST_PDB)
        cls.base_count = surface_voxels(cls.base)
        extreme = cls.base.get_extreme_values()
        cls.far = (extreme[1] + 8, (extreme[2] + extreme[3]) / 2, (extreme[4] + extreme[5]) / 2)

    def added_voxels(self, elem, name):
        atom = self.base.atomList[0].copy()
        atom.elem, atom.atom_name, atom.res, atom.res_no, atom.chain = elem, name, "LIG", 999, "Z"
        atom.record_name = "HETATM"
        atom.x, atom.y, atom.z = self.far
        structure = PDBTools.BioPy_Structure(list(self.base.atomList) + [atom])
        return surface_voxels(structure) - self.base_count

    def check(self, elem, name, voxel_radius):
        self.assertEqual(self.added_voxels(elem, name), self.SPHERE_VOXELS[voxel_radius],
                         "%s/%s" % (elem, name))

    def test_protein_elements_unchanged(self):
        self.check("C", "C1", 3)
        self.check("N", "N1", 2)
        self.check("O", "O1", 2)
        self.check("S", "SD", 2)

    def test_ligand_elements(self):
        self.check("P", "P", 3)
        self.check("SE", "SE", 3)
        self.check("F", "F1", 2)
        self.check("CL", "CL", 3)

    def test_metals_get_default(self):
        self.check("ZN", "ZN", 3)
        self.check("CA", "CA", 3)
        self.check("NA", "NA", 3)

    def test_hydrogens_skipped(self):
        self.check("H", "H1", 0)
        self.check("D", "D1", 0)

    def test_name_fallback_without_element(self):
        self.check("", "CA", 3)
        self.check("", "P", 3)
        self.check("", "OG", 2)

    def test_atom_element(self):
        atom = self.base.atomList[0].copy()
        atom.elem, atom.atom_name = "ZN", "ZN"
        self.assertEqual(GridTools.atom_element(atom), "ZN")
        atom.elem = ""
        self.assertEqual(GridTools.atom_element(atom), "Z")
        atom.elem, atom.atom_name = " se ", "SE"
        self.assertEqual(GridTools.atom_element(atom), "SE")


class LigandOnPath(unittest.TestCase):
    """Examples/test_ligand.pdb: a phosphate on the path between LYS 26 and LYS 125."""

    @classmethod
    def setUpClass(cls):
        cls.without, cls.origin_without = lysine_sasds(LIGAND_PDB, ligand=False)
        cls.with_ligand, cls.origin_with = lysine_sasds(LIGAND_PDB, ligand=True)

    def test_phosphate_is_the_only_hetero_residue(self):
        structure = PDBTools.read_PDB_file(LIGAND_PDB, hetatm=True)
        hetero = [a for a in structure.atomList if a.record_name == "HETATM"]
        self.assertEqual([(a.res, a.res_no, a.atom_name) for a in hetero],
                         [("PO4", 900, n) for n in ("P", "O1", "O2", "O3", "O4")])
        self.assertEqual(len(structure), 1006)

    def test_without_ligand_matches_test_pdb(self):
        reference, _ = lysine_sasds(TEST_PDB, ligand=False)
        self.assertEqual(self.without, reference)

    def test_grid_does_not_move(self):
        self.assertEqual(self.origin_without, self.origin_with)

    def test_ligand_only_lengthens(self):
        self.assertEqual(set(self.without), set(self.with_ligand))
        changed = {k for k in self.without if abs(self.without[k] - self.with_ligand[k]) > 1e-9}
        shorter = {k for k in changed if self.with_ligand[k] < self.without[k]}
        self.assertEqual(shorter, set())
        self.assertGreaterEqual(len(changed), 5)

    def test_pair_26_125(self):
        self.assertGreater(self.with_ligand[(26, 125)] - self.without[(26, 125)], 2.0)


class CommandLine(unittest.TestCase):

    def test_ligand_flag(self):
        plain_dir, _ = run_jwalk(LIGAND_PDB)
        self.addCleanup(shutil.rmtree, plain_dir)
        ligand_dir, _ = run_jwalk(LIGAND_PDB, "-ligand")
        self.addCleanup(shutil.rmtree, ligand_dir)
        plain = read_crosslink_list(plain_dir, LIGAND_PDB)
        ligand = read_crosslink_list(ligand_dir, LIGAND_PDB)
        self.assertEqual(len(plain), 66)
        self.assertEqual(set(plain), set(ligand))
        pair = ("LYS-26-A-CA", "LYS-125-A-CA")
        self.assertGreater(ligand[pair] - plain[pair], 2.0)
        self.assertFalse(any(ligand[k] < plain[k] - 1e-9 for k in plain))

    def test_help_lists_ligand(self):
        result = subprocess.run([sys.executable, JWALK_SCRIPT, "-h"],
                                env=dict(os.environ, PYTHONPATH=SRC_DIR), capture_output=True, text=True)
        self.assertIn("-ligand", result.stdout)


class RSAParser(unittest.TestCase):

    def test_freesasa_2_0_layout(self):
        line = "RES LYS A  21    81.56  42.4  80.74  55.0   0.82   1.8  49.51  42.9  32.06  41.5"
        self.assertEqual(SurfaceTools.parse_rsa_line(line), (21, "A", "LYS", 42.4))

    def test_freesasa_2_1_layout(self):
        line = "RES LYS   A  21    81.56  40.7  80.74  50.9   0.82   2.0  49.51  42.7  32.06  38.0"
        self.assertEqual(SurfaceTools.parse_rsa_line(line), (21, "A", "LYS", 40.7))

    def test_blank_chain(self):
        line = "RES LYS      21    81.56  40.7  80.74  50.9   0.82   2.0  49.51  42.7  32.06  38.0"
        self.assertEqual(SurfaceTools.parse_rsa_line(line), (21, " ", "LYS", 40.7))

    def test_hetero_residue_has_no_relative_sasa(self):
        line = "RES SE4   A 147    37.34   N/A  37.34   N/A   0.00   N/A   0.00   N/A   0.00   N/A"
        self.assertIsNone(SurfaceTools.parse_rsa_line(line))

    def test_insertion_code_and_negative_number(self):
        line = "RES LYS   A  52A   81.56  40.7  80.74  50.9   0.82   2.0  49.51  42.7  32.06  38.0"
        self.assertEqual(SurfaceTools.parse_rsa_line(line)[0], 52)
        line = "RES LYS   A  -3    81.56  40.7  80.74  50.9   0.82   2.0  49.51  42.7  32.06  38.0"
        self.assertEqual(SurfaceTools.parse_rsa_line(line)[0], -3)

    def test_other_lines(self):
        self.assertIsNone(SurfaceTools.parse_rsa_line("REM RES _ NUM      All-atoms   Total-Side"))
        self.assertIsNone(SurfaceTools.parse_rsa_line("RES LYS A 21"))
        self.assertIsNone(SurfaceTools.parse_rsa_line(""))

    def test_bundled_example(self):
        """Examples/test.rsa is FreeSASA 2.0 output: 93 of 124 residues above 7%, all 12 lysines."""
        exposed = set()
        with open(os.path.join(EXAMPLES, "test.rsa")) as rsa:
            for line in rsa:
                parsed = SurfaceTools.parse_rsa_line(line)
                if parsed and parsed[3] > 7.0:
                    exposed.add(parsed[:3])
        self.assertEqual(len(exposed), 93)
        self.assertEqual(sum(1 for r in exposed if r[2] == "LYS"), 12)

    def test_write_pdb_without_water(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        out = os.path.join(directory, "nowater.pdb")
        SurfaceTools.write_pdb_without_water(TEST_PDB, out)
        with open(out) as f:
            lines = f.read().splitlines()
        self.assertFalse(any(l.startswith("HETATM") and l[17:20] == "HOH" for l in lines))
        self.assertEqual(sum(1 for l in lines if l.startswith("HETATM")), 18)
        # 1002 ATOM lines in the file; read_PDB_file drops one hydrogen to give 1001 atoms
        self.assertEqual(sum(1 for l in lines if l.startswith("ATOM")), 1002)


@unittest.skipUnless(shutil.which("freesasa"), "freesasa not on PATH")
class Freesasa(unittest.TestCase):

    def test_surface_runs(self):
        directory, stdout = run_jwalk(TEST_PDB, "-surface")
        self.addCleanup(shutil.rmtree, directory)
        self.assertIn("12 are on the surface", stdout)
        self.assertEqual(len(read_crosslink_list(directory, TEST_PDB)), 66)

    def test_surface_with_ligand_sees_ligand_not_water(self):
        directory, stdout = run_jwalk(LIGAND_PDB, "-surface", "-ligand")
        self.addCleanup(shutil.rmtree, directory)
        with open(os.path.join(directory, "test_ligand.rsa")) as f:
            rsa = f.read()
        self.assertIn("test_ligand_nowater.pdb", rsa)
        self.assertIn("RES PO4", rsa)
        self.assertNotIn("RES HOH", rsa)
        self.assertFalse(os.path.exists(os.path.join(directory, "test_ligand_nowater.pdb")))
        self.assertIn("12 are on the surface", stdout)

    def test_stale_rsa_is_replaced(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        shutil.copy(TEST_PDB, directory)
        with open(os.path.join(directory, "test.rsa"), "w") as f:
            f.write("REM  stale\n")
        env = dict(os.environ, PYTHONPATH=SRC_DIR)
        subprocess.run([sys.executable, JWALK_SCRIPT, "-i", "test.pdb", "-surface", "-ncpus", "1"],
                       cwd=directory, env=env, capture_output=True, text=True, check=True)
        with open(os.path.join(directory, "test.rsa")) as f:
            self.assertNotIn("stale", f.readline())


if __name__ == "__main__":
    unittest.main()
