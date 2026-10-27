# Ligand mode (`-ligand`)

How Jwalk treats ligands, ions and other `HETATM` records, the `-ligand` flag that turns this on, how each atom is turned into an obstacle on the grid, what was wrong before, how the changes were validated, and what is still left undone. Written 2026-10-27 on the `optimize` branch.

## Background

Jwalk reads a PDB file with Biopython and separates the residues into three groups: protein (`ATOM` records), hetero residues (`HETATM` records: ligands, ions, cofactors, modified residues such as MSE) and waters. `PDBTools.read_PDB_file(filename, hetatm=False, water=False)` has always had the switches to include the last two groups, but `jwalk` never passed them, so hetero atoms were always dropped. A ligand sitting in a channel between two lysines was therefore invisible to the path search and the SASD could run straight through it.

Users who wanted ligands counted had been editing `PDBTools.py` by hand to change the `hetatm=False` default. That is what `-ligand` replaces.

## Usage

```
jwalk -i model.pdb -ligand
jwalk -i model.pdb -ligand -xl_list crosslinks.txt
jwalk -i model.pdb -ligand -aa1 LYS -aa2 HIS
```

`-ligand` combines with every other option. Without it, behaviour is exactly as before (verified byte-for-byte, see [Validation](#validation)).

### What is included

With `-ligand`, every residue whose Biopython hetfield starts with `H` is read: anything written as `HETATM` that is not water. That covers ligands, ions, cofactors, buffer molecules (for example the BME and selenate in `Examples/test.pdb`) and modified amino acids (MSE, LLP, SEP, ...).

### What is still excluded

* **Water** (hetfield `W`, residue names HOH, WAT, DOD). Including crystallographic waters as obstacles would change results a great deal and is not what ligand mode is for. `read_PDB_file` has a separate `water=` switch; no command line flag exposes it.
* **Hydrogens and deuteriums** on ligands. Protein hydrogens were already dropped by name in `read_PDB_file`; ligand hydrogens are read but given no sphere (see below), so neither contributes to the surface.
* Anything after the first `CONECT` or `END` record in the file: Biopython stops reading coordinates there. Ligands placed after the `CONECT` block (unusual, but possible in hand-edited files) are silently ignored.

### What `-ligand` does *not* change

* Crosslink sites are still chosen by residue name against `-aa1`/`-aa2` or the `-xl_list`. A modified lysine written as `HETATM` (for example LLP) becomes an obstacle but is not treated as a LYS for crosslinking.
* `get_euclidean_distances` re-reads the file and takes only `ATOM` records named `CA`, so a calcium ion (`HETATM`, name `CA`) is never mistaken for an alpha carbon.
* The `-surface` option runs freesasa on the original file without `--hetatm`, so freesasa's solvent-accessibility test still ignores ligands even when the grid includes them. See [Not done](#not-done).

## How atoms become obstacles

`GridTools.generate_solvent_accessible_surface` places a voxel sphere around every atom (except the side chains of the residues being crosslinked). The radius of the sphere is the atom's radius plus 0.8 A (a probe allowance Jwalk has always used), divided by the voxel size and rounded to a whole number of voxels. Any voxel touched by at least one sphere is inside the protein and the breadth-first search may not enter it.

### Before this change

The sphere was chosen from the **first letter of the atom name**: `C` -> carbon sphere, `N` -> nitrogen, `O` -> oxygen, `S` -> sulfur, anything else -> no sphere at all. That is right for every atom of the twenty standard amino acids, whose names always start with the element letter, and it was all Jwalk needed while hetero atoms were never read. With ligands it goes wrong:

* **Ignored entirely:** P, Br, I, and every metal (Zn, Fe, Mg, Mn, K, Cu, Ni, Co, ...). A phosphate, a nucleotide's phosphates, or a zinc ion added nothing to the surface.
* **Wrong sphere:** Se (`SE`) got the sulfur sphere, Na (`NA`) got the nitrogen sphere, Cl (`CL`) and Ca (`CA`) got the carbon sphere (which happens to be the right size at 1 A voxels, by luck).

There was a second, hidden problem. `BioPyAtom` tried to record the element with `atom.get_element()`, but that method does not exist in Biopython (the element is the `.element` attribute). The call raised `AttributeError`, a bare `except` swallowed it, and `elem` was `""` for every atom. So even if the sphere code had looked at the element, there was nothing to look at.

### After this change

1. `BioPyAtom.elem` is now the element Biopython read from columns 77-78 of the PDB record, upper-cased. When that column is blank Biopython infers the element from the atom name using the PDB alignment convention, which correctly tells `CA ` (calcium, left-justified in columns 13-14) from ` CA` (alpha carbon, in columns 14-15).

2. `GridTools.atom_element(atom)` returns that element, or, if it is still empty, the first letter of the atom name (the old rule), so files with neither an element column nor an inferable name behave exactly as before.

3. The sphere is chosen by element from this table:

   | element | radius (A) | radius + 0.8 | voxels at `-vox 1` | voxels at `-vox 2` | source |
   |---|---|---|---|---|---|
   | C | 1.73 | 2.53 | 3 | 1 | unchanged Jwalk value |
   | N | 1.43 | 2.23 | 2 | 1 | unchanged Jwalk value |
   | O | 1.30 | 2.10 | 2 | 1 | unchanged Jwalk value |
   | S | 1.67 | 2.47 | 2 | 1 | unchanged Jwalk value |
   | P | 1.80 | 2.60 | 3 | 1 | Bondi 1964 |
   | Se | 1.90 | 2.70 | 3 | 1 | Bondi 1964 |
   | F | 1.47 | 2.27 | 2 | 1 | Bondi 1964 |
   | Cl | 1.75 | 2.55 | 3 | 1 | Bondi 1964 |
   | Br | 1.85 | 2.65 | 3 | 1 | Bondi 1964 |
   | I | 1.98 | 2.78 | 3 | 1 | Bondi 1964 |
   | H, D | - | - | none | none | not placed |
   | any other element | 1.73 | 2.53 | 3 | 1 | **default: the carbon radius** |

   The default is deliberately the largest of the protein radii. An atom Jwalk does not know (a metal ion, an unusual element) gets a sphere at least as big as a carbon's, so it is never left out of the surface again. That slightly over-blocks small ions (a bare Zn2+ is smaller than a carbon) but errs on the side of a blocked path rather than a path through an atom. If a specific element needs a different radius, add it to the `radius` dictionary in `generate_solvent_accessible_surface`.

   Note the rounding: at `-vox 1` there are only two sphere sizes in practice (2 or 3 voxels, 33 or 123 voxels in volume), and at `-vox 2` every atom gets the same 1-voxel sphere. The Bondi values are on a slightly different scale from the four historical Jwalk values (Bondi S is 1.80, Jwalk's is 1.67), which is why P rounds to 3 voxels while S rounds to 2. This only matters at `-vox 1`.

4. Spheres are built once per distinct voxel radius and the expansion loop is unchanged, so the speed of the surface step is the same as before.

The four protein elements keep their original radii, and every standard amino-acid atom name starts with its element letter, so **for protein atoms the new rule picks the same sphere as the old one**. That is why runs without `-ligand` are unaffected.

## Effects on results

### Ligands only ever lengthen true shortest paths

Adding obstacles can only remove routes, so the true shortest path between two residues can only stay the same or get longer. Verified with a Dijkstra search on the same grids: with `-ligand` on `Examples/test.pdb`, no pair's true shortest path got shorter.

### But Jwalk's SASD can occasionally get shorter

Jwalk's search is a first-arrival breadth-first search: it keeps the path with the fewest voxel steps and reports its length in angstroms. That is not the same as the shortest path in angstroms, and it is not monotonic. Blocking a voxel can remove a hop-minimal route that used many diagonal steps and replace it with a route one step longer that is more axial and shorter in angstroms. On `test.pdb` with the grid held fixed, 2 of 66 lysine pairs got shorter with `-ligand` (by 0.5 A and 0.6 A) for exactly this reason. This is the original Jwalk algorithm, not something introduced here; it is written up in `20261026-performance-results.md` under "SASDs are hop-minimal, not angstrom-minimal".

### The grid moves

`GridTools.makeGrid` sizes the grid from the extreme coordinates of every atom read. A ligand that sticks out past the protein enlarges the grid and shifts its origin, which changes how every atom falls into voxels. On `test.pdb` the selenate shifts the origin by 0.29 A in y, and that alone changed 64 of 66 SASDs between a plain run and a `-ligand` run, by amounts of the same order as the voxel size. Do not compare individual SASDs between a run with `-ligand` and one without and attribute the difference to the ligand; most of it is discretisation. Sizing the grid from protein atoms only would remove this effect and would be a small change to `makeGrid`, but it has not been made.

### Numbers on the bundled example

`Examples/test.pdb` (1FGA, basic fibroblast growth factor) has two selenate ions (SE4 147 and 301), two beta-mercaptoethanols (BME 302 and 303) and 65 waters. `-ligand` reads 18 extra atoms (1019 instead of 1001). Relative to the first `-ligand` implementation (name-based spheres), the element-based spheres change 4 of 66 SASDs, all longer by 0.1-1.4 A, all involving residues 26, 119, 125 or 21 next to SE4 301: the two Se atoms moved from the 2-voxel sulfur sphere to a 3-voxel one. Nothing else in the file is affected because C, O and S atoms map the same way under both rules.

## Validation

All runs used `Jwalk/src/jwalk-launcher` (Python 3.12, NumPy 2.5.4, Biopython 1.88) on `Examples/test.pdb` unless stated. Outputs were compared with `cmp`/`md5sum` on both `Jwalk_results/*_crosslink_list.txt` and `*_crosslinks.pdb`.

**Determinism.** Four `-ligand` runs (three with the default CPU count, one with `-ncpus 1`) gave identical output, so hash comparisons are meaningful.

**No change without `-ligand`.** Against the code before each change, these runs were byte-identical in both output files: default (`-vox 1`), `-vox 2`, `-xl_list Examples/xl_list.txt`, `-aa1 HIS -aa2 LYS`. A copy of `test.pdb` with every `HETATM` line removed gave identical output with and without `-ligand`.

**Atoms read.** `read_PDB_file(..., hetatm=True)` returns 1019 atoms versus 1001, the extra 18 being SE4 x2 (5 atoms each) and BME x2 (4 atoms each). No waters.

**Per-element spheres.** A synthetic atom of each element was appended to the structure in empty space 8 A beyond the protein, and the number of surface voxels it added was compared with the expected sphere volume (123 voxels for a 3-voxel radius, 33 for 2, 0 for none). All 21 cases matched:

| case | old voxels | new voxels |
|---|---|---|
| C, N, O, S (protein elements) | 123, 33, 33, 33 | same |
| Se | 33 | 123 |
| P, Br, I | 0 | 123 |
| F | 0 | 33 |
| Cl, Ca | 123 | 123 |
| Zn, Fe, Mg, K | 0 | 123 |
| Na | 33 | 123 |
| H, D | 0 | 0 |
| no element, name `CA` | 123 | 123 (fallback: carbon) |
| no element, name `P` | 0 | 123 |
| no element, name `NA` | 33 | 33 (fallback: nitrogen, unchanged limitation) |

**End to end with an atom that used to be invisible.** A single `HETATM ... P PO4` atom was placed on a voxel of the 26-125 path from a `-ligand` run (inserted before the `CONECT` block, after the last `HETATM`). With the new code, `-ligand` changed 9 of 66 SASDs, all longer; 26-125 went from 25.91 A to 28.86 A. With the previous code (`git stash`) the same file gave output identical to the run without the phosphorus: the atom was ignored.

**Monotonicity.** On fixed grids (ligand atoms added to the surface, grid origin held), Dijkstra shortest paths for all 66 pairs never got shorter with ligands; 13 changed. Jwalk's BFS changed 13 and shortened 2 (see above). Jwalk's BFS overestimates Dijkstra by 1.1 A on average, 6.3 A at most, and never underestimates, with or without ligands. These measurements were taken with the first, name-based `-ligand` implementation; the element change alters the Se spheres and so would move the exact figures slightly, but not the conclusions.

**Not tested.** `-surface` (freesasa is not installed on this machine). `setup.py install`. Files whose element column is absent *and* whose atom names Biopython cannot interpret.

## Not done

* **freesasa ignores ligands under `-surface`.** `SurfaceTools.check_solvent_accessibility_freesasa` runs `freesasa --format=rsa --radii=naccess file.pdb` without `--hetatm`. With `-ligand -surface`, residue accessibility is judged without the ligand while the path search includes it. Adding `--hetatm` when `-ligand` is set is a one-line change, but freesasa's naccess radius set may not cover every ligand element and this needs testing on a machine with freesasa.
* **`SurfaceTools.check_solvent_accessibility`** (the pure-Python accessibility check, not the freesasa one) still indexes radii by the first letter of the atom name and would raise `KeyError` on a phosphorus atom. Nothing in `jwalk` calls it.
* **Grid sized from all atoms.** See "The grid moves" above.
* **No `-water` flag.** `read_PDB_file(water=True)` exists; nothing exposes it.
* **Radii for metals.** All use the carbon default. A per-element table (ionic or vdW) could be added to the `radius` dictionary if someone needs it.
* **No ligand-in-a-channel example file.** The phosphorus test above was built by hand; `Examples/` still has no file where a ligand sits on a crosslink path.

## Reproducing

```
cd Jwalk
python3.12 -m venv src/venv && src/venv/bin/pip install -r requirements.txt
mkdir -p /tmp/jw && cp Examples/test.pdb Examples/xl_list.txt /tmp/jw && cd /tmp/jw
/path/to/Jwalk/src/jwalk-launcher -i test.pdb            # 66 SASDs
/path/to/Jwalk/src/jwalk-launcher -i test.pdb -ligand    # 66 SASDs, 4 differ from the name-based version
```

To check that an atom was actually read before interpreting a result:

```python
from Jwalk import PDBTools
s = PDBTools.read_PDB_file('test.pdb', hetatm=True)
print([(a.res, a.res_no, a.atom_name, a.elem) for a in s.atomList if a.record_name == 'HETATM'])
```
