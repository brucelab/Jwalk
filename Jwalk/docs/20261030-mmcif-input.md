# mmCIF input

Jwalk reads mmCIF files (`.cif`, `.mmcif`) directly, next to PDB files. Written 2026-10-30 on the
`optimize` branch after reviewing two xlinkdb structures, `7ohr.cif` (51 chains) and `9syg.cif`
(84 chains, two-character chain ids), with the aim of no longer converting mmCIF to PDB first.

## Usage

```
jwalk -lys -i 9syg.cif
jwalk -lys -i 9syg.cif -ligand
jwalk                              # every .pdb, .cif and .mmcif file in the directory
```

Every option works as for a PDB file. Results go to `Jwalk_results/<base name>_crosslink_list.txt`
and `<base name>_crosslinks.pdb`, so `9syg.cif` gives `9syg_crosslink_list.txt` with `9syg.cif` in
the Model column. Chain ids are written as they are in the file (`LYS-27-A0-CA`); the `MODEL` line of
the path file runs residue and chain together (`LYS27A0-LYS36A0`) as it always has, which is
ambiguous with chain ids longer than one character, so parse the crosslink list rather than that.

## Why

xlinkdb converts mmCIF files to PDB before running Jwalk. The PDB format allows 62 one-character
chain ids, so entries with more chains (9syg has 84) are skipped, and the conversion renames the
chains, which has to be undone afterwards with a `_chainnames.txt` file. The review also found that
the conversion loses atoms: it renames every nucleotide to PHE, its phosphorus to CA, and strips
the primes from the sugar atom names (C2', C4', C5'), which then collide with the base atoms C2, C4
and C5 of the same residue, and Biopython's PDB parser drops the duplicates. On 7OHR the three RNA
chains lose 15% of their atoms, 9,604 heavy atoms, which leaves holes in the rRNA that the path
search goes through.

7OHR, lysine pairs, 6 CPUs, same grid origin on both routes (chains mapped back with the
`_chainnames.txt` file):

| | PDB route (converted file) | mmCIF route |
|---|---|---|
| SASDs within 60 A | 50,893 | 42,416 |
| pairs found on one route only | 8,567 | 90 |
| common pairs with a different SASD | 10,420 of 42,326 (25%), median 0.8 A | |
| of which longer on the mmCIF route | 8,353, up to 55 A (LYS 112 - LYS 118 of chain C: 11.7 -> 66.7 A) | |
| wall time | 34.0 s | 34.5 s |

Every production Jwalk result for an RNA-containing complex made through the conversion carries
this error.

## What is read

Reading is done with Biopython's `MMCIFParser`, which gives the structure the author identifiers
(`auth_asym_id`, `auth_seq_id`, the ones a PDB file carries), so chains and residue numbers are
the same as in the equivalent PDB and nothing downstream of `read_PDB_file` changed. Parsing
7ohr.cif (236k atom rows) takes 4.9 s against 2.0 s for the PDB; the whole run is unchanged because
the parsed dictionary is reused for the polymer lookup below.

* **Polymer chains, always.** Every residue of a polymer entity is read, including the modified
  residues that an mmCIF file writes as `HETATM` records: MSE, SEP, TPO, and in RNA PSU, OMG, A2M,
  OMC, OMU, 5MC, ... The file says which residues belong to a chain (`_pdbx_poly_seq_scheme`,
  matched on author chain, residue number and insertion code), so these are told from ligands.
  They are obstacles like any other residue. They are not crosslink sites: a modified lysine
  written as MLY or LLP is not a LYS, exactly as before.
* **Ligands and ions only with `-ligand`**, as for PDB files. Water never.
* **Hydrogens** are dropped by element (H, D). PDB files are still judged by the first letter of the
  atom name, as before, so their output is unchanged byte for byte.
* **Alternative locations**: Biopython keeps one, as for PDB files. Only the first model is read.
* **Insertion codes** are still ignored in residue keys, as everywhere in Jwalk.

A PDB file has no record of which `HETATM` residues belong to a chain (Biopython does not expose
`MODRES`), so for PDB input modified residues still count as ligands. That is the one difference in
meaning between the two formats.

### Why polymer HETATM residues matter

9syg has 3,138 atoms in modified rRNA nucleotides. Left out, they are holes in the rRNA; counted
as ligands they come and go with the ions. Lysine pairs at 6 CPUs:

| 9syg | SASDs | vs the row above |
|---|---|---|
| modified nucleotides dropped (as a plain hetero-residue rule would) | 42,458 | |
| modified nucleotides included (now the default) | 42,007 | 465 pairs gone, 14 new; 2,059 changed, 1,705 longer, up to 48 A |
| plus `-ligand` (197 Mg, 80 K, 24 Na, 6 Zn) | 41,623 | 396 gone, 12 new; 1,759 changed, 1,430 longer, up to 30 A |

Run time 40 s and 1.1 GB per worker for each.

## Euclidean distances

`get_euclidean_distances` used to re-read the input file with fixed PDB columns. It now takes the
C-alpha coordinates from the structure already in memory (any format), rounded to the three
decimals of the file so that PDB output is byte-identical to before. The old file reading is kept
for callers that pass no structure.

## `-surface` with mmCIF

freesasa reads mmCIF with `--cif`. Jwalk gives it a copy of the file (`<name>_jwalk.cif`, removed
afterwards) holding exactly the atoms of its grid: the chains with their modified residues,
ligands only with `-ligand`, no water, and runs it with `--hetatm` so that the `HETATM` rows that
are left are counted. For PDB input nothing changed.

In a `.rsa` file made from mmCIF, freesasa writes the chain label and residue number as one word
(`A022` is residue 22 of chain `A0`), which cannot be split without knowing the chains. The residues
of interest are therefore looked up by that word (`SurfaceTools.rsa_words`), first as written and
then without a trailing insertion code. The PDB layouts of FreeSASA 2.0 and 2.1 parse as before.

## Validation

* `tests/test_mmcif.py`, 23 tests (47 with `tests/test_ligand.py`), about 6 s, run with
  `src/venv/bin/python -m unittest discover -s tests -v` from the `Jwalk` directory. They use
  `Examples/test.cif` and `Examples/test_modres.cif`, made from `Examples/test.pdb` by
  `tests/make_test_cif.py` with the full RCSB column set and a `_pdbx_poly_seq_scheme`;
  `test_modres.cif` writes MET 142 as MSE.
* `test.cif` gives the same 1,001 atoms as `test.pdb`, the same 18 ligand atoms with `-ligand`,
  and the `jwalk` output differs only in the file name in the Model column; the path file is
  identical. `test_modres.cif` reads the 8 MSE atoms without `-ligand`, with element SE.
* On the production structures (6 CPUs): the new code on `7OHR.pdb` is byte-identical to the
  previous commit; `7ohr.cif` runs in 34.5 s (42,416 SASDs); `9syg.cif` in 40.4 s (42,007 SASDs,
  74 chains with lysines, two-character ids throughout).
* `-surface` on `test.cif` gives the same 66 SASDs and "12 are on the surface" as on `test.pdb`;
  with `-ligand` the `.rsa` lists SE4 and BME and no water; on `test_modres.cif` it lists MSE.
  freesasa 2.1.3 reads both the RCSB files and the generated examples; it warns about ligand
  atoms whose element it guesses, as before.

## Not done

* **Output naming.** Results take the file's base name as given, so `9syg.cif` gives lower-case
  `9syg_...`; xlinkdb's result reader and its chain-name mapping step expect the renamed chains of
  the converted PDB and would need to take the real chain ids instead.
* **Grid sized from every atom read**, as before; see the ligand-mode notes.
* **Metal ions under `-ligand`** still get the carbon radius. On xlinkdb-written PDB files the
  element column is absent and iron is inferred as fluorine from the name; mmCIF files carry the
  element, so this does not arise with mmCIF input.
* **`SurfaceTools.check_solvent_accessibility`** (the pure-Python check nothing calls) is untouched.
