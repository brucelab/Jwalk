# Jwalk performance results

Results of implementing the changes proposed in `20261026-performance.md`, one step at a time on the `optimize` branch. After every step the program was run on a set of test cases and its output compared with the previous step; the results below are from one final back-to-back measurement of the original code and every step.

## Summary

* **Large benchmark** (`big.pdb`, test.pdb tiled into 4 chains, 4,008 atoms, 356 SASDs, 4 CPUs): **33.7 s -> 0.73 s wall (46x)** and **116 s -> 3.1 s CPU (38x)**.
* **`-xl_list` with 124 crosslinks from one residue**: 30.2 s -> 0.53 s wall (57x).
* **Bundled example** (`Examples/test.pdb`, `-ncpus 4`): 1.57 s -> 0.31 s wall. The remaining time is mostly importing NumPy and Biopython (about 1.4 s of CPU), which no step here touches.
* Output is byte-identical between consecutive steps on every test case, with **one deliberate exception**: step 1 fixes a bug in the original path extraction (see [Behaviour change](#behaviour-change-a-bug-fix-in-step-1)).

## Step by step

Times are the best of 3 runs (the baseline benchmark is a single run because it is slow). "Step gain" is the CPU-time ratio against the previous step.

| Step | Commit | Change | Benchmark wall | Benchmark CPU | Step gain (CPU) | Total (CPU) |
|---|---|---|---:|---:|---:|---:|
| 0 | original | Python 3 port of the original code | 33.73 s | 116.5 s | | 1.0x |
| 1 | `4751170` | BFS rewritten with flat arrays, parent pointers, precomputed step costs (items 1, 9) | 4.46 s | 15.1 s | 7.71x | 7.7x |
| 2 | `d3feb15` | Grid and voxel data sent to workers once via Pool initializer; pool closed (item 3) | 4.16 s | 14.9 s | 1.01x | 7.8x |
| 3 | `a03adc2` | `xl_list`: one BFS per start residue; BFS stops once all targets are reached (item 2) | 4.17 s | 15.0 s | 1.00x | 7.8x |
| 4 | `380ab2f` | Level-synchronous NumPy BFS | 1.45 s | 4.7 s | 3.21x | 24.9x |
| 5 | `09981fc` | First-path-to-a-voxel found with `np.minimum.at` instead of `np.unique` | 1.08 s | 3.4 s | 1.39x | 34.7x |
| 6 | `4de176a` | Euclidean distances only for the SASD pairs, not all N x N (item 7) | 0.90 s | 3.2 s | 1.06x | 36.9x |
| 7 | `b1f201f` | Sphere stamping in `generate_solvent_accessible_surface` vectorized (item 5) | 0.71 s | 3.0 s | 1.04x | 38.2x |
| 8 | `4a1efe5` | Set lookups instead of list scans; search shells built once (item 6) | 0.73 s | 3.1 s | 0.98x | 37.6x |
| 9 | `8f1adc5` | No duplicate atom construction on read; direct formatting of path PDB records (item 8) | 0.73 s | 3.1 s | 1.01x | 37.9x |

## All test cases, original vs final

Wall / CPU seconds, best of 3 (original: best of 3 except the benchmark).

| Case | What it runs | Original wall | Final wall | Original CPU | Final CPU |
|---|---|---:|---:|---:|---:|
| `lys` | all lysine pairs, 4 CPUs | 1.57 | 0.31 | 6.50 | 1.66 |
| `lys_1cpu` | all lysine pairs, 1 CPU | 4.79 | 0.33 | 6.02 | 1.57 |
| `xl` | `-xl_list`, 3 crosslinks, 1 CPU | 1.33 | 0.23 | 2.54 | 1.47 |
| `xl_4cpu` | `-xl_list`, 3 crosslinks, 4 CPUs | 0.72 | 0.29 | 2.74 | 1.55 |
| `xl_shared` | `-xl_list`, 33 lines sharing start residues, duplicates and reversed pairs, 4 CPUs | 3.87 | 0.30 | 14.43 | 1.60 |
| `xl_shared_1cpu` | same list, 1 CPU | 11.78 | 0.28 | 13.02 | 1.50 |
| `xl_errors` | `-xl_list`, 128 lines incl. missing residue/chain, 2 CPUs | 30.16 | 0.53 | 59.97 | 2.02 |
| `asp_lys` | `-aa1 ASP -aa2 LYS` | 1.03 | 0.29 | 3.51 | 1.57 |
| `cys_glu_d30` | `-aa1 CYS -aa2 GLU -max_dist 30` | 0.58 | 0.27 | 1.97 | 1.52 |
| `lys_vox2` | `-vox 2` | 0.43 | 0.28 | 2.08 | 1.56 |
| `big_lys` | large benchmark, 4 CPUs | 33.73 | 0.73 | 116.47 | 3.07 |

## Behaviour change: a bug fix in step 1

The original `calculate_SASDs` / `calculate_specific_SASD` extracted the shortest path by modifying the shared path list and distance of each end voxel in place (`insert(0, ...)`, `append(...)`, `distance[...] +=`). When two end residues share a voxel, the second pair therefore got a path with a duplicated C-alpha at each end and an inflated distance, and could choose a longer route. Step 1 computes both locally, so each SASD is the true shortest path.

How this was verified: a test-only copy of the new BFS that re-creates the old in-place mutation produced output byte-identical to the original on every test case, so the only difference between the original and step 1 is this bug.

Effect on the test data (original -> step 1):

| Case | SASDs | Distances that changed | Paths that changed |
|---|---:|---:|---:|
| `lys` (`test.pdb`) | 66 | 1 (LYS-125 / LYS-135: 16.226 A -> 15.948 A) | 2 |
| `cys_glu_d30` | 15 | 1 | 1 |
| `lys_vox2` | 66 | 3 | 7 |
| `big_lys` | 356 | 4 | 10 |
| `asp_lys`, all `-xl_list` cases | | 0 | 1 in `asp_lys`, 0 in the rest |

Changed distances were always shorter than the original value. The duplicated C-alpha entries in `*_crosslinks.pdb` for those pairs are also gone. Results from the original code for any structure with crosslinks sharing end voxels may therefore be slightly overestimated.

## How correctness was checked

* A harness ran 11 cases per step: all lysine pairs (1 and 4 CPUs), `-xl_list` (1 and 4 CPUs), an `-xl_list` with shared starts, duplicates and reversed pairs, an `-xl_list` with missing residues and chains, other residue types, `-max_dist 30`, `-vox 2`, plus the large benchmark. For each case it compared the stdout, the SASD list (distances to 1e-6) and every path in the `_crosslinks.pdb` file with the reference.
* Steps 2-9 were additionally compared **byte for byte** with the previous step. Step 1 is the reference; see above for why it differs from the original.
* Steps 6-9 touch code that the large cases only partly exercise, so they were also checked directly against the previous version: the solvent accessible grid is exactly equal element for element (step 7); `mark_CAlphas_pairs` on three crosslink lists including bad residues, `find_surface_voxels` at five voxel sizes with and without `surface`, and 200 random `update_crosslink_pairs` cases give identical return values and identical printed text (step 8); the written path file is identical (step 9).

## Fixes made after the optimization steps

Checking the finished work for problems turned up one regression introduced by the optimization (memory) and three crashes in the original code. Each was fixed in its own commit, after the nine steps above, and each was checked by re-running all test cases (byte-identical output to the previous commit) plus a test aimed at the fix. The timings above were not re-measured after these commits; one run of the large benchmark afterwards was within the noise of the step 9 figures.

| Commit | Problem | Fix | Specific test |
|---|---|---|---|
| `6e30983` | **Memory (a regression from step 4).** The NumPy BFS held about 25 bytes per voxel over the whole padded grid, for every search in every worker: 640 MB per worker at 300^3 voxels and 1.4 GB at 400^3, and its `int32` indices would overflow beyond about 1290^3. | A search can't get further than `max_dist / vox` voxels (plus the step that finds a voxel) from its start voxels, so only that window of the grid is held. End voxels outside the window can't have been reached. Indices switch to `int64` if a window ever exceeded `int32`. Peak memory for one search from the centre: 300^3 640 MB -> 131 MB, 400^3 1370 MB -> 164 MB. | 4 new cases with a small `-max_dist` (window smaller than the grid); 600 random searches compared voxel by voxel (visited, distance, parent, root) with the previous implementation, 444 with a window smaller than the grid. The old search never reached a voxel outside the new window. |
| `6222180` | **`-i dir/file.pdb` crashed** with `FileNotFoundError` after the whole calculation, because output names were built with `pdb[:-4]`, which keeps the directory. | Name the outputs from the base name without extension. They still go in `./Jwalk_results`. | Relative and absolute input paths give the same results as running in the file's directory. |
| `64f8e06` | **A blank line in the `-xl_list` file** was reported as a formatting error and aborted the run. | Skip blank and whitespace-only lines; the line counter still counts them, so real errors report the right line. | A list with blank, space-only and tab-only lines gives byte-identical results to the clean list; a malformed line still stops at the right line number. |
| `a7cf7ec` | **`-vox 0.5` crashed** with a `ValueError` traceback and `-vox 0` would divide by zero. | `-vox`, `-max_dist` and `-ncpus` are `type=int` in argparse, so a bad value gives a usage error (exit status 2); a voxel size below 1 is rejected. | Six bad values rejected, two valid ones still run. |

Two behaviours to be aware of: because outputs are named from the file's base name, two inputs with the same file name in different directories now write to the same result files; and voxel size is still whole numbers only (decimal sizes were not tested through the rest of the code, so they are rejected rather than accepted).

### Dead code removed in step 8

In `SurfaceTools.update_crosslink_pairs` a `buried_residues` list was appended to but never read. One of its checks (`x2 in remove_aa2`) compared a 3-tuple with a list of 2-lists and could never be true. It was removed with the set conversion in commit `4a1efe5`; the function was compared with the old version on 200 random inputs, including the text it prints.

### Known issues not fixed

Verified by running them:

* A non-numeric residue number in the `-xl_list` file (for example `21|A|oops|A|`) still stops with a raw `ValueError` traceback instead of the formatting message.
* `PDBTools.calculate_centre_of_mass` raises `UnboundLocalError` for a structure with no atoms (the divisions are also inside the loop). Nothing in `jwalk` uses the result.
* `GridTools.mapGridPosition` returns `0` for an atom outside the grid and `markCAlphas` then fails with `TypeError`. The grid is built around the structure, so this was only reproduced by moving the grid origin.

Found by reading the code only (not run, since freesasa is not installed here, or not tested):

* The `-surface` message prints `len(aa1_CA)-1` lysines "and 1 N-terminus", but nothing in the code adds an N-terminus, so the lysine count shown is one too low.
* `check_solvent_accessibility_freesasa` has `test_open.close` without parentheses (the file is not explicitly closed) and `sd_res` is unbound if the residue list is empty.
* Residues are keyed by number and chain, ignoring insertion codes, in `GridTools` and `get_euclidean_distances`; `get_euclidean_distances` reads every model of an NMR file while Biopython reads only the first.

### SASDs are hop-minimal, not angstrom-minimal (original algorithm, not changed)

The BFS in `SASDTools` is a first-arrival search: the path to a voxel is the one with the fewest steps, and the SASD is the angstrom length (1, sqrt 2, sqrt 3 per step) summed along that path. This is not a weighted shortest path. Measured against Dijkstra on the same 26-neighbour grid (`Examples/test.pdb`, 66 LYS pairs, `-vox 1`, 2026-10-27): the BFS overestimates by 1.1 A on average and 6.3 A at most, and never underestimates.

A consequence is that adding obstacles can make an SASD *shorter*. With `-ligand` on `test.pdb` and the grid held fixed, 13 of 66 SASDs changed and 2 got shorter (26-125: 24.02 -> 23.56 A; 46-77: 56.35 -> 55.71 A). In both, a blocked hop-minimal route was replaced by a route one hop longer but with fewer diagonal steps. Dijkstra on the same grids never got shorter. This reproduces the original queue BFS, so it is not an optimization regression; replacing it with Dijkstra would shift every result by about 1 A and is left as a separate decision.

Note also that `-ligand` can move the grid origin (`makeGrid` sizes the grid from every atom read), which changes the voxel discretisation for every residue. On `test.pdb` that alone changed 64 of 66 SASDs; hold the grid fixed when comparing runs with and without ligands.

These measurements were taken with the first `-ligand` implementation, which chose each atom's sphere from the first letter of its name. The sphere is now chosen by element (see `20261027-ligand-mode.md`), which changes the two Se spheres on `test.pdb` and would move the exact figures slightly, but not the conclusions.

### Fixed since: element never read from the PDB

`BioPyAtom` recorded the element with `atom.get_element()`, which does not exist in Biopython (the element is the `.element` attribute); a bare `except` hid the `AttributeError` and `elem` was `""` for every atom. Nothing in `jwalk` used the element until the surface code was switched to element-based radii, so this had no visible effect before. Fixed in `PDBTools.BioPyAtom.__init__`; the element column of PDBs written by `BioPy_Structure.write_to_PDB` is now filled in too (`write_sasd_to_pdb` does not use that writer, so result files are unchanged).

## Caveats

* **Shared machine.** The server had a load average of about 10-12 during the final measurement, and earlier in the session the same code varied by up to 1.5x between runs (the original benchmark ranged from 34 s to 78 s wall). CPU time is steadier than wall time but not immune. Treat the ratios as indicative; the ordering of the steps is reliable, the exact factors are not.
* **Test structures are small.** All inputs are `Examples/test.pdb` (1,002 atoms, one chain) or four translated copies of it. Larger real complexes will spend a bigger share of time in the BFS and may scale differently.
* `-surface` (needs freesasa) and `setup.py install` were not run.
* Small steps (7, 8, 9) are lost in the noise of whole-program timings; their individual effect was measured on the function (step 7: `generate_solvent_accessible_surface` 311 ms -> 69 ms; step 9: `read_PDB_file` 37.7 -> 29.0 ms and `write_sasd_to_pdb` 40.6 -> 31.6 ms, on the large structure).

## Proposals not implemented

* **Item 4, `SurfaceTools.check_solvent_accessibility` (cKDTree).** The `jwalk` command never calls it; `-surface` uses the freesasa variant. It is dead code on the command-line path, so there was nothing to measure, and its float32/float64 arithmetic would make an exact-equivalence check fragile.
* `calculate_centre_of_mass` and `get_extreme_values` in `PDBTools.py`: not hot spots, and vectorizing the centre-of-mass sum would change the order of floating point additions.
* **Compiled BFS (Numba/Cython).** Not needed after step 4-5, and it would add a build dependency.

## Where the time goes now

A profile of the large benchmark with 1 CPU after step 9 puts about 65% of the run in `breadth_first_search`; most of the rest is importing NumPy and Biopython, parsing the PDB and writing results. Any further speedup would have to come from the BFS itself, for example by sharing work between searches or compiling it.

## Reproducing

The harness (test cases, runner, comparison scripts) was kept in the session scratch directory and is not part of this repository. Each step is a separate commit on `optimize`, so any step can be checked out and run against `Examples/test.pdb` with the commands in the README.
