# Jwalk performance review

Reviewed files: `src/Jwalk/` (`SASDTools.py`, `GridTools.py`, `SurfaceTools.py`, `PDBTools.py`) and the `jwalk` launcher. `venv/`, `build/` and `py2_original/` were skipped. Nothing was profiled; the ranking comes from reading the code and the complexity of each loop. Profile on a real PDB before committing to a fix order.

## High impact

**1. The BFS in `SASDTools.py` is the main cost.**
It runs from every start residue, up to `max_dist`, over the whole grid, in `calculate_SASDs` (L179) and `calculate_specific_SASD` (L29).
- **`queue.pop(0)` (L113, L261):** on a list this is O(n) per pop. Use `collections.deque` with `popleft()`.
- **Full path copy per voxel (L122, L273):** `temp_list = visited[x_n,y_n,z_n][:]` copies the whole path for every voxel reached. That makes time and memory quadratic in path length. Store a parent pointer per voxel instead, and rebuild the path only for the few end voxels needed.
- **Tuple-keyed dicts:** `visited` and `distance` are dicts keyed by `(x,y,z)` tuples. Replace them with flat NumPy arrays indexed by `z*ny*nx + y*nx + x`. Cache `x_size()`/`y_size()`/`z_size()` in locals, because they are method calls inside the innermost loop.
- **Distance branching (L281-293):** the if/elif chain per neighbour is redundant. Precompute a step cost for each of the 26 entries in `comb`, e.g. `comb = [(dx, dy, dz, cost), ...]`.
- **Best option:** write the BFS in Numba or Cython, or use `scipy.sparse.csgraph.dijkstra` on a graph of empty voxels. Either should give one to two orders of magnitude.

**2. The BFS does repeated work in the `xl_list` mode.**
`calculate_specific_SASD` runs a full BFS for each crosslink pair. Group the pairs by start residue and run one BFS per unique start. Also stop early once all target end voxels have been reached, instead of going all the way to `max_dist`.

**3. `parallel_BFS` serializes the grid for every task (L386-422).**
- **Repeated pickling:** `itertools.repeat(dens_map)`, `aa1_voxels` and the other arguments are pickled with each task chunk, which repeatedly copies the full 3D array.
- **Fix:** use a `Pool(initializer=...)` that stores these as module globals once per worker, or rely on fork inheritance. Pass only the start residue to `map`.
- **Cleanup:** the pool is never closed or joined, so use `with Pool(...)`.
- **Scheduling:** consider `imap_unordered` with a `chunksize`, since BFS cost varies a lot per residue.

**4. `check_solvent_accessibility` in `SurfaceTools.py` (L~95) is O(residue atoms × 30 points × all atoms) in pure Python.**
The code's own comment says "this is not very efficient". Each sphere point is tested against every atom in the structure.
- **Fix:** put all atom coordinates in a NumPy array and build a `scipy.spatial.cKDTree`. Query only the neighbours within the maximum radius of each point (about 3.5 Å). That makes it a near-linear lookup.
- **Smaller issue:** it also loops over every atom of each residue of interest, not just CA. Check whether that is intended.

## Medium impact

**5. `generate_solvent_accessible_surface` in `GridTools.py` (L~233).**
The sphere is stamped one voxel at a time in Python, with a bounds check per voxel and chained `[z][y][x]` indexing.
- Precompute each sphere as an `(n,3)` offset array.
- Clip to the grid and add with NumPy slicing (`fullMap[z0:z1, y0:y1, x0:x1] += sphere_mask[...]`), or use `np.add.at`.
- The four near-identical `if/elif` blocks (C/O/N/S) collapse into one dict lookup.
- The `mapGridPosition()` call per atom can be vectorised across all atoms at once.

**6. List membership tests in `GridTools.py` and `SurfaceTools.py`.**
- **`mark_CAlphas_pairs`:** `[atom.res_no, atom.chain] in aa1` is a linear list scan for every atom, so cost is atoms × crosslinks. Use sets of tuples. The same applies to `x not in atom_check` and `i not in index_to_delete`.
- **`find_surface_voxels`:** the sphere shells are built with `[x,y,z] not in sphereN` list scans. Build each shell by distance bands instead: `r_prev² < d² <= r²`. It is computed twice per run (once for aa1, once for aa2), so hoist or cache it.
- **`SurfaceTools.update_crosslink_pairs`:** the same `i not in index_to_delete` and `x not in buried_residues` list pattern. Use sets.

**7. `get_euclidean_distances` (`SASDTools.py` L440) builds an N×N dict.**
It computes the Euclidean distance between every pair of CA atoms. That is millions of entries for large complexes, but only the pairs in `sasds` are looked up. Look up the two CA coordinates per SASD directly. It also re-parses the PDB file with string slicing, when the structure is already in memory.

## Lower impact

**8. `PDBTools.py`**
- **`read_PDB_file`:** every atom is built twice. `BioPyAtom(atom)` is called as a bare statement and then again for `append`. Remove the first call. Biopython's `PDBParser` is also slow, so `Bio.PDB.FastMMCIFParser` or a minimal fixed-width parser would help.
- **`calculate_centre_of_mass`:** the three `x_CoM/y_CoM/z_CoM` divisions are inside the loop. They should be after it. It is wasteful, and it also gives the wrong result for an empty list.
- **`get_extreme_values`:** this makes about 12 getter calls per atom. Use `min()`/`max()` over a coordinate array.
- **Atom storage:** `numpy.append` on object arrays copies the whole array. A coordinate/mass NumPy array stored on the structure would speed up many of the steps above.
- **`write_sasd_to_pdb`:** it builds a full `BioPyAtom` object per path voxel just to format one line. Format the string directly and write in one batch with `writelines`.

**9. A path-extraction oddity in `SASDTools.py` (L154-163, L309-318).**
- **What it does:** `visited[(x,y,z)].insert(0, ...)` and `.append(...)` mutate the shared path lists, and `distance[(x,y,z)] += ...` mutates shared state, once per end residue.
- **Why it matters:** that is a possible correctness bug if two end residues share a voxel. It also means `for i in [1, len(...)-1]` only adds two segments; this is probably intentional but unclear.
- **What to do:** avoid mutating `visited` and `distance`. Compute the final distance locally. Parent-pointer reconstruction (item 1) removes the problem.

## Suggested order
1. Do the BFS rewrite (items 1-2) and the Pool initializer (item 3). Together they will dominate runtime.
2. Add the cKDTree solvent accessibility (item 4) if you use the non-`-surface` path on large structures.
3. Do the set and NumPy clean-ups (items 5-8). They are easy wins.

There are no tests in the repo, so before changing anything, save the output of `Examples/test.pdb` as a baseline and confirm results stay identical afterwards.
