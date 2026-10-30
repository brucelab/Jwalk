#===============================================================================
#     This file is part of Jwalk.
#
#     Tests for mmCIF input: reading .cif files, modified residues of a chain
#     (HETATM records that belong to a polymer) always included, ligands only
#     with -ligand, the Euclidean distances taken from the structure, results
#     identical to the PDB version of the same structure, long chain ids, and
#     -surface (freesasa --cif) on mmCIF input.
#
#     Examples/test.cif and Examples/test_modres.cif are made from
#     Examples/test.pdb by tests/make_test_cif.py.
#
#     Run from the Jwalk directory with the launcher's virtual environment:
#
#         src/venv/bin/python -m unittest discover -s tests -v
#
#     The -surface tests are skipped when freesasa is not on the PATH.
#===============================================================================

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
TEST_CIF = os.path.join(EXAMPLES, "test.cif")
MODRES_CIF = os.path.join(EXAMPLES, "test_modres.cif")
JWALK_SCRIPT = os.path.join(SRC_DIR, "Jwalk", "jwalk")

sys.path.insert(0, SRC_DIR)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Jwalk import PDBTools, GridTools, SurfaceTools, SASDTools
from test_ligand import run_jwalk, read_crosslink_list, lysine_sasds


def atom_keys(structure):
    return [(a.record_name, a.chain, a.res_no, a.res, a.atom_name) for a in structure.atomList]


def result_text(directory, name):
    with open(os.path.join(directory, "Jwalk_results", name + "_crosslink_list.txt")) as f:
        return f.read()


class ReadMMCIF(unittest.TestCase):

    def test_is_mmcif(self):
        self.assertTrue(PDBTools.is_mmcif("x.cif"))
        self.assertTrue(PDBTools.is_mmcif("dir/9SYG.CIF"))
        self.assertTrue(PDBTools.is_mmcif("x.mmcif"))
        self.assertFalse(PDBTools.is_mmcif("x.pdb"))
        self.assertFalse(PDBTools.is_mmcif("x.cif.pdb"))

    def test_same_atoms_as_pdb(self):
        pdb = PDBTools.read_PDB_file(TEST_PDB)
        cif = PDBTools.read_PDB_file(TEST_CIF)
        self.assertEqual(len(pdb.atomList), 1001)
        self.assertEqual(atom_keys(pdb), atom_keys(cif))
        for a, b in zip(pdb.atomList, cif.atomList):
            self.assertAlmostEqual(a.x, b.x, places=3)
            self.assertAlmostEqual(a.y, b.y, places=3)
            self.assertAlmostEqual(a.z, b.z, places=3)
            self.assertEqual(a.elem, b.elem)

    def test_ligands_only_with_hetatm_and_never_water(self):
        without = PDBTools.read_PDB_file(TEST_CIF)
        with_ligands = PDBTools.read_PDB_file(TEST_CIF, hetatm=True)
        self.assertEqual(sum(1 for a in without.atomList if a.record_name == "HETATM"), 0)
        self.assertEqual(len(with_ligands.atomList) - len(without.atomList), 18)
        names = set(a.res for a in with_ligands.atomList if a.record_name == "HETATM")
        self.assertEqual(names, {"SE4", "BME"})
        self.assertEqual(atom_keys(with_ligands), atom_keys(PDBTools.read_PDB_file(TEST_PDB, hetatm=True)))

    def test_polymer_residues(self):
        residues = PDBTools.polymer_residues(MMCIF2Dict(TEST_CIF))
        self.assertEqual(len(residues), 124)
        self.assertIn(("A", 21, " "), residues)
        self.assertIn(("A", 142, " "), residues)
        self.assertNotIn(("A", 147, " "), residues)  # SE4 147 is a ligand
        self.assertEqual(PDBTools.polymer_residues({}), set())

    def test_modified_residue_always_read(self):
        plain = PDBTools.read_PDB_file(TEST_PDB)
        modres = PDBTools.read_PDB_file(MODRES_CIF)
        self.assertEqual(len(modres.atomList), len(plain.atomList))
        mse = [a for a in modres.atomList if a.res_no == 142]
        self.assertEqual(len(mse), 8)
        self.assertTrue(all(a.res == "MSE" and a.record_name == "HETATM" for a in mse))
        self.assertEqual([a.elem for a in mse if a.atom_name == "SE"], ["SE"])
        # still no ligands without hetatm
        self.assertEqual(set(a.res for a in modres.atomList if a.record_name == "HETATM"), {"MSE"})
        # and with hetatm the ligands come too
        self.assertEqual(len(PDBTools.read_PDB_file(MODRES_CIF, hetatm=True).atomList), len(plain.atomList) + 18)

    def test_modified_residue_is_an_obstacle_not_a_site(self):
        structure = PDBTools.read_PDB_file(MODRES_CIF)
        grid = GridTools.makeGrid(structure, 1)
        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "LYS", "LYS")
        self.assertEqual(len(aa1_CA), 12)
        self.assertNotIn((142, "A", "MSE"), aa1_CA)
        aa1_CA, aa2_CA = GridTools.markCAlphas(grid, structure, "MET", "LYS")
        self.assertNotIn((142, "A", "MSE"), aa1_CA)  # written as MSE, so not a MET either

    def test_hydrogens_dropped_by_element(self):
        # a hydrogen whose name starts with a digit is kept by the PDB rule but dropped by element
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        path = os.path.join(directory, "h.cif")
        with open(TEST_CIF) as f:
            text = f.read()
        extra = ("ATOM 9999 H 1HB . LYS A 1 2 ? 27.0 36.0 14.0 1.00 10.00 ? 21 LYS A 1HB 1\n"
                 "ATOM 9998 C CX . LYS A 1 2 ? 27.0 36.0 15.0 1.00 10.00 ? 21 LYS A CX 1\n")
        text = text.replace("ATOM   1 ", extra + "ATOM   1 ", 1)
        with open(path, "w") as f:
            f.write(text)
        names = [a.atom_name for a in PDBTools.read_PDB_file(path).atomList if a.res_no == 21]
        self.assertIn("CX", names)
        self.assertNotIn("1HB", names)


class WriteSubset(unittest.TestCase):

    def test_atoms_only(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        out = os.path.join(directory, "subset.cif")
        PDBTools.write_mmcif_subset(TEST_CIF, out, lambda group, chain, number, icode, name: group == "ATOM")
        with open(out) as f:
            lines = f.readlines()
        with open(TEST_CIF) as f:
            original = f.readlines()
        self.assertEqual(sum(1 for l in lines if l.startswith("HETATM")), 0)
        # every ATOM row is kept (1002 rows: one atom has two alternative locations)
        self.assertEqual(sum(1 for l in lines if l.startswith("ATOM")), sum(1 for l in original if l.startswith("ATOM")))
        self.assertEqual([l for l in lines if not l.startswith(("ATOM", "HETATM"))],
                         [l for l in original if not l.startswith(("ATOM", "HETATM"))])
        self.assertEqual(len(PDBTools.read_PDB_file(out, hetatm=True).atomList), 1001)

    def test_for_surface(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        out = os.path.join(directory, "surface.cif")
        PDBTools.write_mmcif_for_surface(MODRES_CIF, out, ligands=False)
        with open(out) as f:
            het = [l.split()[17] for l in f if l.startswith("HETATM")]
        self.assertEqual(set(het), {"MSE"})
        PDBTools.write_mmcif_for_surface(MODRES_CIF, out, ligands=True)
        with open(out) as f:
            het = [l.split()[17] for l in f if l.startswith("HETATM")]
        self.assertEqual(set(het), {"MSE", "SE4", "BME"})


class Euclidean(unittest.TestCase):

    def test_structure_matches_file(self):
        structure = PDBTools.read_PDB_file(TEST_PDB)
        sasds = {((21, "A", "LYS"), (26, "A", "LYS"), 10.0): [], ((26, "A", "LYS"), (125, "A", "LYS"), 24.0): []}
        from_file = SASDTools.get_euclidean_distances(sasds, TEST_PDB, "LYS", "LYS")
        from_structure = SASDTools.get_euclidean_distances(sasds, TEST_PDB, "LYS", "LYS", structure)
        self.assertEqual(from_file, from_structure)
        from_cif = SASDTools.get_euclidean_distances(sasds, TEST_CIF, "LYS", "LYS", PDBTools.read_PDB_file(TEST_CIF))
        self.assertEqual(from_file, from_cif)
        self.assertEqual(len(from_file), 2)

    def test_unknown_residue_dropped(self):
        structure = PDBTools.read_PDB_file(TEST_PDB)
        sasds = {((21, "A", "LYS"), (999, "A", "LYS"), 10.0): []}
        self.assertEqual(SASDTools.get_euclidean_distances(sasds, TEST_PDB, "LYS", "LYS", structure), {})


class CommandLine(unittest.TestCase):

    def test_cif_gives_pdb_results(self):
        pdb_dir, _ = run_jwalk(TEST_PDB)
        self.addCleanup(shutil.rmtree, pdb_dir)
        cif_dir, stdout = run_jwalk(TEST_CIF)
        self.addCleanup(shutil.rmtree, cif_dir)
        self.assertIn("66 SASDs calculated", stdout)
        # the files differ only in the input file name in the Model column
        self.assertEqual(result_text(pdb_dir, "test"), result_text(cif_dir, "test").replace("test.cif     ", "test.pdb     "))
        with open(os.path.join(pdb_dir, "Jwalk_results", "test_crosslinks.pdb")) as a, \
             open(os.path.join(cif_dir, "Jwalk_results", "test_crosslinks.pdb")) as b:
            self.assertEqual(a.read(), b.read())

    def test_modified_residue_runs(self):
        directory, stdout = run_jwalk(MODRES_CIF)
        self.addCleanup(shutil.rmtree, directory)
        self.assertIn("66 SASDs calculated", stdout)
        self.assertEqual(len(read_crosslink_list(directory, MODRES_CIF)), 66)

    def test_default_file_list_includes_cif(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        shutil.copy(TEST_CIF, directory)
        env = dict(os.environ, PYTHONPATH=SRC_DIR)
        result = subprocess.run([sys.executable, JWALK_SCRIPT, "-ncpus", "1"], cwd=directory, env=env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("calculating crosslinks on test.cif", result.stdout)
        self.assertTrue(os.path.exists(os.path.join(directory, "Jwalk_results", "test_crosslink_list.txt")))

    def test_help_mentions_mmcif(self):
        result = subprocess.run([sys.executable, JWALK_SCRIPT, "-h"], env=dict(os.environ, PYTHONPATH=SRC_DIR),
                                capture_output=True, text=True)
        self.assertIn("mmCIF", result.stdout)


class LongChainIds(unittest.TestCase):

    def test_output_with_two_character_chains(self):
        directory = tempfile.mkdtemp(prefix="jwalk_test_")
        self.addCleanup(shutil.rmtree, directory)
        cwd = os.getcwd()
        os.chdir(directory)
        try:
            path = [[1, 1, 1], [2, 2, 2], [3, 3, 3]]
            sasds = {((27, "A0", "LYS"), (36, "Az", "LYS"), 50.0, 28.5): path}
            PDBTools.write_sasd_to_txt(sasds, "9syg.cif")
            grid = GridTools.Map(__import__("numpy").zeros((5, 5, 5)), [0, 0, 0], 1, "g")
            PDBTools.write_sasd_to_pdb(grid, sasds, "9syg.cif")
        finally:
            os.chdir(cwd)
        text = result_text(directory, "9syg")
        self.assertIn("LYS-27-A0-CA", text)
        self.assertIn("LYS-36-Az-CA", text)
        with open(os.path.join(directory, "Jwalk_results", "9syg_crosslinks.pdb")) as f:
            self.assertTrue(f.readline().startswith("MODEL 1 LYS27A0-LYS36Az"))

    def test_rsa_words(self):
        residues = {(22, "A0", "GLU"): None, (148, "Az", "LYS"): None, (217, "B", "ALA"): None}
        words = SurfaceTools.rsa_words(residues)
        self.assertEqual(words[("GLU", "A022")], (22, "A0", "GLU"))
        self.assertEqual(SurfaceTools.parse_rsa_line("RES GLU  A022    111.28  64.7  45.52  34.9  65.75 157.6  60.15 100.3  51.13  45.6", words),
                         (22, "A0", "GLU", 64.7))
        self.assertEqual(SurfaceTools.parse_rsa_line("RES LYS  Az148   241.12 120.4 169.60 107.0  71.52 161.4 169.60 110.2  71.52 155.2", words),
                         (148, "Az", "LYS", 120.4))
        # insertion code on the residue number
        self.assertEqual(SurfaceTools.parse_rsa_line("RES ALA   B217A    0.05   0.0   0.00   0.0   0.05   0.1   0.00   0.0   0.05   0.1", words),
                         (217, "B", "ALA", 0.0))
        # a residue that is not of interest, and a hetero residue without relative SASA
        self.assertIsNone(SurfaceTools.parse_rsa_line("RES GLY  A023    20.72  19.2  11.27  17.6   9.45  21.5  11.30  15.9   9.42  25.7", words))
        self.assertIsNone(SurfaceTools.parse_rsa_line("RES PSU  AA61    164.53   N/A  83.09   N/A  81.44   N/A  21.0   N/A  10.0   N/A", words))
        # the PDB layouts still parse when words is given
        self.assertEqual(SurfaceTools.parse_rsa_line("RES LYS A  21   102.55  42.4  85.10  54.0  17.45  17.0  58.87  73.7  43.68  26.9", words),
                         (21, "A", "LYS", 42.4))


@unittest.skipUnless(shutil.which("freesasa"), "freesasa is not on the PATH")
class SurfaceMMCIF(unittest.TestCase):

    def test_surface_on_cif_matches_pdb(self):
        pdb_dir, pdb_out = run_jwalk(TEST_PDB, "-surface")
        self.addCleanup(shutil.rmtree, pdb_dir)
        cif_dir, cif_out = run_jwalk(TEST_CIF, "-surface")
        self.addCleanup(shutil.rmtree, cif_dir)
        self.assertIn("12 are on the surface", cif_out)
        self.assertEqual(read_crosslink_list(pdb_dir, TEST_PDB), read_crosslink_list(cif_dir, TEST_CIF))
        with open(os.path.join(cif_dir, "test.rsa")) as f:
            rsa = f.read()
        self.assertIn("RES LYS", rsa)
        self.assertNotIn("SE4", rsa)
        self.assertNotIn("HOH", rsa)
        self.assertFalse(os.path.exists(os.path.join(cif_dir, "test_jwalk.cif")))

    def test_surface_with_ligand_on_cif(self):
        directory, stdout = run_jwalk(TEST_CIF, "-surface", "-ligand")
        self.addCleanup(shutil.rmtree, directory)
        with open(os.path.join(directory, "test.rsa")) as f:
            rsa = f.read()
        self.assertIn("SE4", rsa)
        self.assertIn("BME", rsa)
        self.assertNotIn("HOH", rsa)
        self.assertIn("12 are on the surface", stdout)

    def test_surface_sees_modified_residue(self):
        directory, stdout = run_jwalk(MODRES_CIF, "-surface")
        self.addCleanup(shutil.rmtree, directory)
        with open(os.path.join(directory, "test_modres.rsa")) as f:
            rsa = f.read()
        self.assertIn("MSE", rsa)
        self.assertNotIn("SE4", rsa)
        self.assertIn("12 are on the surface", stdout)


if __name__ == "__main__":
    unittest.main()
