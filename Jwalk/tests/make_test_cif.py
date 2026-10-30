"""
Writes the mmCIF versions of Examples/test.pdb used by the tests:

    Examples/test.cif         the same atoms as test.pdb, in mmCIF form
    Examples/test_modres.cif  the same, with MET 142 written as the modified residue MSE
                              (HETATM records, SD -> SE), as deposited structures do

Both carry the _pdbx_poly_seq_scheme category that marks the residues of the polymer chain, and
the full set of _atom_site columns (label_* and auth_*) that RCSB files have, which freesasa
needs to read the file with --cif.

Run from the Jwalk directory:  src/venv/bin/python tests/make_test_cif.py
"""
import os

JWALK_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PDB = os.path.join(JWALK_DIR, "Examples", "test.pdb")

COLUMNS = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
           "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
           "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_formal_charge",
           "auth_seq_id", "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]

SCHEME = ["asym_id", "entity_id", "seq_id", "mon_id", "ndb_seq_num", "pdb_seq_num",
          "auth_seq_num", "pdb_mon_id", "auth_mon_id", "pdb_strand_id", "pdb_ins_code", "hetero"]


def quote(name):
    return '"%s"' % name if "'" in name else name


def read_pdb(path):
    atoms = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            atoms.append(dict(
                group=line[:6].strip(), serial=int(line[6:11]), name=line[12:16].strip(),
                altloc=line[16].strip() or ".", resname=line[17:20].strip(), chain=line[21],
                resseq=int(line[22:26]), icode=line[26].strip() or "?",
                x=line[30:38].strip(), y=line[38:46].strip(), z=line[46:54].strip(),
                occ=line[54:60].strip(), b=line[60:66].strip(),
                element=(line[76:78].strip() or line[12:16].strip()[0]).upper()))
    return atoms


def write_cif(atoms, path, title):
    # polymer residues of chain A in order, label_seq_id 1..n; each hetero residue is its own
    # label_asym_id and entity, as in RCSB files
    polymer = []
    for a in atoms:
        if a["group"] == "ATOM" or a["modres"]:
            key = (a["chain"], a["resseq"], a["icode"])
            if key not in [p[0] for p in polymer]:
                polymer.append((key, a["resname"]))
    seq_id = dict((key, i + 1) for i, (key, name) in enumerate(polymer))
    # each hetero residue gets its own label_asym_id (B, C, ... Z, AA, AB, ...) and one entity per
    # hetero residue type, numbered from 2, as in RCSB files
    letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    labels = list(letters[1:]) + [x + y for x in letters for y in letters]
    het_asym = {}
    het_types = sorted(set(a["resname"] for a in atoms if a["group"] == "HETATM" and not a["modres"]))
    rows = []
    for a in atoms:
        key = (a["chain"], a["resseq"], a["icode"])
        if key in seq_id:
            asym, entity, label_seq = "A", 1, str(seq_id[key])
        else:
            if key not in het_asym:
                het_asym[key] = labels[len(het_asym)]
            asym = het_asym[key]
            entity = 2 + het_types.index(a["resname"])
            label_seq = "."
        rows.append([a["group"], str(a["serial"]), a["element"], quote(a["name"]), a["altloc"],
                     a["resname"], asym, str(entity), label_seq, a["icode"], a["x"], a["y"], a["z"],
                     a["occ"], a["b"], "?", str(a["resseq"]), a["resname"], a["chain"],
                     quote(a["name"]), "1"])
    widths = [max(len(r[i]) for r in rows) for i in range(len(COLUMNS))]
    with open(path, "w") as f:
        f.write("data_test\n#\n_entry.id   test\n#\n_struct.title   '%s'\n#\nloop_\n" % title)
        for c in SCHEME:
            f.write("_pdbx_poly_seq_scheme.%s\n" % c)
        for (key, name) in polymer:
            f.write("A 1 %d %s %d %d %d %s %s %s . n\n" % (seq_id[key], name, seq_id[key], key[1], key[1], name, name, key[0]))
        f.write("#\nloop_\n")
        for c in COLUMNS:
            f.write("_atom_site.%s\n" % c)
        for r in rows:
            f.write(" ".join(v.ljust(w) for v, w in zip(r, widths)).rstrip() + "\n")
        f.write("#\n")


def main():
    atoms = read_pdb(PDB)
    for a in atoms:
        a["modres"] = False
    write_cif(atoms, os.path.join(JWALK_DIR, "Examples", "test.cif"),
              "Examples/test.pdb (1FGA) as mmCIF")

    # MET 142 -> MSE: selenomethionine, written as HETATM with SE in place of SD
    for a in atoms:
        if a["group"] == "ATOM" and a["chain"] == "A" and a["resseq"] == 142:
            a["group"] = "HETATM"
            a["modres"] = True
            a["resname"] = "MSE"
            if a["name"] == "SD":
                a["name"] = "SE"
                a["element"] = "SE"
    write_cif(atoms, os.path.join(JWALK_DIR, "Examples", "test_modres.cif"),
              "Examples/test.pdb (1FGA) as mmCIF with MET 142 as the modified residue MSE")


if __name__ == "__main__":
    main()
