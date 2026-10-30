# Jwalk benchmark: original vs optimized (2026-10-30)

Structures: xlinkdb's `Writeonly_jwalk/complLigMemory_102426_1/` set of PDB files, those whose names
begin with "1" (11YB, 11AA, 12AH, 11HV, 11HE), plus 7OHR for both codes at 6 CPUs. The same directory
holds the production Jwalk results of Oct 26-27 referred to below (`Jwalk_results/` and
`Ligand/Jwalk_results/`). Mode: `jwalk -lys -i <pdb> -ncpus 16` (same as production, which runs plain
`jwalk -lys`, apart from the CPU count). Machine: the xlinkdb server, 24 cores, idle apart from these
runs. One interpreter for every variant: the venv of the production install (Python 3.12.15,
NumPy 2.5.4, Biopython 1.88).

## Variants

| Variant | Sources | Notes |
|---|---|---|
| `orig` | the production install of the original (the `jwalk` on PATH) | Python 3 port of the original. Identical to brucelab commit `ba84b51` after reverting one local edit (see caveat 1). |
| `opt_bugs` | brucelab `2c4d84e` + bug emulation | All speed optimizations, with the original path-extraction bug re-created (`patch_opt_bugs.py`, see Files). Nothing else differs from `opt_fixed`. |
| `opt_fixed` | brucelab `2c4d84e` (HEAD) | Optimized, bugs fixed. |

The first speed commit in brucelab (`4751170`) also contains the bug fix, so no commit is "speed only";
`opt_bugs` is a copy of HEAD whose `calculate_SASDs` mutates the shared end-voxel path and distance in
place exactly as the original did.

This report measures the finished code only. The step-by-step history of the nine optimization commits,
the other fixes made afterwards and how each was verified, and the design notes on the hop-minimal
search, ligands and the grid origin are in `20261026-performance-results.md`; `-ligand` itself is
described in `20261027-ligand-mode.md`.

## Run time (16 CPUs; optimized = best of 2 runs, original = 1 run)

| File | Atoms | SASDs | orig wall | opt_bugs wall | opt_fixed wall | Speedup | orig CPU | opt_fixed CPU | CPU speedup |
|---|---|---|---|---|---|---|---|---|---|
| 11YB | 25,509 | 2,290 | 68.0 s | 3.1 s | 2.9 s | 23x | 958 s | 34 s | 28x |
| 11AA | 29,520 | 11,366 | 171.3 s | 8.8 s | 8.1 s | 21x | 2,332 s | 92 s | 25x |
| 12AH | 60,310 | 24,689 | 340.7 s | 18.6 s | 16.3 s | 21x | 4,185 s | 181 s | 23x |
| 11HV | 147,551 | 39,522 | 734.7 s | 42.8 s | 34.9 s | 21x | 9,002 s | 401 s | 22x |
| 11HE | 147,580 | 39,436 | 734.6 s | 41.4 s | 34.1 s | 22x | 8,828 s | 393 s | 22x |
| Total | | | 2,049 s (34 min) | 115 s | 96 s | 21x | 25,305 s | 1,102 s | 23x |

CPU time is user+system over all worker processes. Repeats of the optimized variants agreed within 1 s.
`opt_bugs` is slightly slower than `opt_fixed` because emulating the bug means building a full path list
for every end voxel, as the original did.

Peak memory (largest single process, i.e. roughly one worker):

| File | orig | opt_fixed |
|---|---|---|
| 11YB | 0.7 GB | 0.15 GB |
| 11AA | 1.7 GB | 0.25 GB |
| 12AH | 10.6 GB | 0.48 GB |
| 11HV | 24.2 GB | 0.82 GB |
| 11HE | 24.2 GB | 0.82 GB |

The original stores a complete copy of the path for every voxel it reaches, so its memory grows with the
square of the search radius; 16 such workers on 29IW/29IX (224k atoms) would be at risk on this machine.

## Results

**Speed optimizations alone change nothing.** `opt_bugs` output is byte-identical to `orig` on all five
structures, both `<pdb>_crosslink_list.txt` and `<pdb>_crosslinks.pdb`. `orig` is in turn byte-identical
to the production results of Oct 26 in `Writeonly_jwalk/complLigMemory_102426_1/Jwalk_results/`.

**Bug fix (orig / opt_bugs -> opt_fixed).** Same number of SASDs, same pairs, same Euclidean distances.
A few percent of the SASDs get shorter; none get longer.

| File | SASDs | Changed | % | Mean change | Median | Largest | Paths changed | Cross 35 A | Cross 42 A | Cross 51 A |
|---|---|---|---|---|---|---|---|---|---|---|
| 11YB | 2,290 | 45 | 2.0 | -0.76 A | -0.52 A | -3.21 A | 91 | 1 | 0 | 1 |
| 11AA | 11,366 | 360 | 3.2 | -0.65 A | -0.45 A | -6.00 A | 687 | 4 | 5 | 5 |
| 12AH | 24,689 | 569 | 2.3 | -0.51 A | -0.39 A | -8.44 A | 1,091 | 1 | 7 | 4 |
| 11HV | 39,522 | 1,174 | 3.0 | -0.53 A | -0.38 A | -5.94 A | 2,587 | 4 | 4 | 12 |
| 11HE | 39,436 | 1,247 | 3.2 | -0.57 A | -0.44 A | -5.78 A | 2,623 | 3 | 7 | 12 |
| 7OHR | 50,893 | 2,038 | 4.0 | -0.53 A | -0.39 A | -5.23 A | 4,137 | 14 | 18 | 33 |

"Cross N A" = pairs whose SASD moves from above N A to at or below it.

More paths change than distances: the original also appended two extra C-alpha records to a path it
had already stored for an earlier pair, so some pairs keep their distance but lose two spurious atoms
in the `_crosslinks.pdb` file.

### The bug

In the original `SASDTools.calculate_SASDs`, after the breadth-first search from a start residue, the
path to each end residue is read from the shared `visited` dictionary. The code inserted the start
C-alpha and appended the end C-alpha into that shared list in place, and added the two C-alpha link
lengths to the shared `distance` entry. When two end residues share a surface voxel (common for
neighbouring lysines), the second one received a path with a duplicated C-alpha at each end and an
inflated distance, and could therefore pick a longer route. The optimized code computes the path and
distance locally per pair. This is the only fix that changes `jwalk -lys` output; the other brucelab
fixes (element reading, memory window, `-i dir/file.pdb`, blank `-xl_list` lines, `-vox` validation,
FreeSASA 2.1 parsing) had no effect here, as the byte-identical `opt_bugs` vs `orig` result shows.

### Reproducing the original's output: `-includebugs`

The optimized code has a flag that gives the original's
output without the original's run time: `jwalk -lys -includebugs`. It switches on the original path
extraction (`SASDTools.original_shortest_paths`, the same emulation the `opt_bugs` variant above
used) and the original's surface spheres chosen by the first letter of the atom name
(`GridTools.generate_solvent_accessible_surface(..., legacy_spheres=True)`), so that with `-ligand`
it also matches the hand-edited production tree of Oct 27 (P and metal atoms get no sphere).
Checked on the benchmark outputs: `-includebugs` on 11YB.pdb is byte-identical to the original run,
and `-includebugs -ligand` on 12AH.pdb is byte-identical to the original with `hetatm=True`
(the Oct 27 production Ligand file). Everything else (speed, memory, mmCIF input, the crash fixes)
stays.

## 7OHR (235,845 atoms, 51 chains, 50,893 SASDs), both codes at 6 CPUs

| Run | CPUs | Wall | CPU (all workers) | Largest process RSS |
|---|---|---|---|---|
| `orig`, measured | 6 | 1,374.6 s (22.9 min) | 6,929 s | 24.8 GB |
| `opt_bugs`, measured | 6 | 51.2 s | 243 s | 1.0 GB |
| `opt_fixed`, measured | 6 | 36.4 s | 152 s | 1.0 GB |
| ratio | | 37.8x | 45.6x | 24x |
| original, production Oct 26 (`Jwalk_results/` file timestamps, 7NAC 18:29:57 -> 7OHR 18:53:03) | 24 (default) | ~23 min 06 s | n/a | n/a |
| original, production Oct 27 "Ligand" run (`Ligand/Jwalk_results/`, 7NAC 22:46:20 -> 7OHR 23:09:16) | 24 (default) | ~22 min 56 s | n/a | n/a |

The measured original run is byte-identical to the Oct 26 production file, and the `opt_bugs` run is
byte-identical to the original run (md5 a9c4226d... for the list, dc6f6269... for the path PDB). The bugs-fixed run shortens
2,038 of the 50,893 SASDs (see the results table above); nothing else differs.

The wall-time ratio is larger here than the 21x seen at 16 CPUs because at 6 CPUs the optimized run is
dominated by single-threaded work (parsing the PDB, building the grid, writing 50,893 paths) and only
keeps 4.2 cores busy on average, while the original's search keeps all 6 busy; the CPU-time ratio is the
like-for-like figure. Note also that the original took the same ~23 min with 6 CPUs as production did
with the default 24, so the original does not scale beyond a handful of workers on a structure this
size (each task re-pickles the whole grid, and the workers need ~25 GB each at peak).

## Ligand mode check on 12AH (added later the same day)

Question: do the Oct 27 "Ligand" production results differ from the optimized `-ligand` output? Tested
on 12AH (60,310 ATOM, 545 HETATM: ADP with 18 P, two FeMo clusters with 90 Fe and 6 Mo, homocitrate,
2-oxoglutarate), all at 16 CPUs.

| Comparison | SASDs | Changed | Shorter | Longer | Largest |
|---|---|---|---|---|---|
| 1. production Oct 27 Ligand vs original with `hetatm=True` forced | 24,398 / 24,398 | 0 (byte-identical) | | | |
| 2. original `hetatm=True` vs `opt_bugs -ligand` (element spheres only) | 24,398 / 24,398 | 59 | 16 | 43 | 1.27 A |
| 3. `opt_bugs -ligand` vs `opt_fixed -ligand` (bug fix only) | 24,398 / 24,398 | 573 | 573 | 0 | 8.44 A |
| 4. production Oct 27 Ligand vs `opt_fixed -ligand` (what a user sees) | 24,398 / 24,398 | 627 | 584 | 43 | 8.44 A |
| 5. production Oct 26 no-ligand vs Oct 27 Ligand (ligands in the original) | 24,689 / 24,398 | 409 of the 24,397 common | 45 | 364 | 49.1 A |

Row 1 shows the Oct 27 Ligand files were made with the original code and the `hetatm=True` edit (that
edit was reverted in the production tree at 12:23 on 2026-10-30). Row 2 is the element-based sphere
change on its own: the original gave no sphere to P, Fe or Mo, so 114 of the 545 hetero atoms were
invisible to the search; counting them lengthens 43 SASDs by up to 1.27 A and shortens 16 (the hop-minimal
BFS effect described in `20261027-ligand-mode.md`). Row 3 is the same end-voxel bug fix as in the
non-ligand runs. Row 5 is for scale: including ligands at all is a far larger effect (292 pairs lost,
364 lengthened, up to 49 A) than either code change.

## Caveats

1. `PDBTools.py` in the production install of the original was edited on 2026-10-27 15:04 so that
   `read_PDB_file` defaults to `hetatm=True`. Every other copy of the original (the Python 2 sources, the
   egg, the build directory, xlinkdb's own copy, brucelab `ba84b51`) defaults to `False`, and the Oct 26
   production results match `False`. The benchmark used `False`, which is the genuine original and also
   the optimized default. The edit was reverted at 12:23 on 2026-10-30, so the `jwalk` on PATH again
   excludes HETATM.
2. Five files begin with "1", not four: 12AH was included.
3. The original was run once per structure; the optimized variants twice (best shown). The machine was
   otherwise idle (load < 0.5 before the runs).
4. GNU time's max RSS is the largest single process, not the sum over workers.
5. `-surface` and `-xl_list` modes were not benchmarked, and `-ligand` only for results on 12AH; the bug
   emulation covers the default all-pairs mode only (with or without `-ligand`).

## Files

Harness and small result files, in `benchmark/` next to this document's parent `Jwalk/` directory:

- Scripts: `run_all.sh` (whole benchmark), `run_variant.sh <variant> <ncpus> <rep>` (one variant over
  `pdb_order.txt`), `patch_opt_bugs.py` (turns a copy of the optimized `SASDTools.py` into `opt_bugs`),
  `compare.py <run_A> <run_B>` (SASD, Euclidean and path comparison of two `Jwalk_results/`),
  `summarize_timings.py`.
- Timings: `timings.log` (one line per run: variant, repeat, file, wall, user, sys, max RSS in kB, exit
  status), `timings_7OHR.log`; per-run `*.log` (jwalk stdout) and `*.time` under `runs/`.
- Comparisons: `cmp_orig_vs_bugs.json`, `cmp_bugs_vs_fixed.json`, `cmp_orig_vs_fixed.json`.

Not kept in the repository: the variant source copies (`orig` is the production install with the one
line reverted, the others are commit `2c4d84e` with and without the patch) and the `Jwalk_results/`
output of every run (about 440 MB, regenerable with the scripts).

To repeat: put the PDB files or symlinks to them in a directory, write their names to `pdb_order.txt`,
set `B` at the top of `run_variant.sh`/`run_all.sh` to the directory holding `src/<variant>/Jwalk` and
`PY` to a Python with NumPy and Biopython, and run `NCPUS=16 ./run_all.sh`.
