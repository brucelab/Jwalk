# Treating chosen chain:residue sites as lysine (proposal)

A proposed option that makes Jwalk treat a user-supplied list of `chain:residue` sites as LYS when it reads a structure. Nothing described here is implemented yet. Written 2026-11-13 after reading the code on the `optimize` branch.

## Why

A user maps crosslinked proteins onto PDB chains by sequence search (BLAST), so the structure is often a homolog or variant of the crosslinked protein. At a crosslinked lysine, the aligned residue in the structure may not be lysine (Arg or Gln, for example). Jwalk's default all-lysine scan only starts from residues named LYS, so those sites are skipped.

The user's current workaround is to rename the residue to LYS at those sites when converting the mmCIF to PDB, and run Jwalk on the converted file. Jwalk now reads mmCIF directly, so the conversion is the only reason left to go through PDB. A substitution list would also leave the structure file untouched.

## What the code does today

* **The start point is the CA atom.** `markCAlphas` (`GridTools.py`) selects residues with `atom.res == aa1` (or `aa2`) and takes only the atom named `CA`. No lysine side-chain atom (NZ) is needed, so renaming the residue to LYS is all that is required.
* **`-xl_list` already works for any residue.** `mark_CAlphas_pairs` (`GridTools.py`) matches on `(res_no, chain)` only and never checks the residue type. A crosslink list can name a non-lysine site today. The problem exists only for the default all-lysine scan (`-lys`, or `-aa1`/`-aa2`).
* **Side chains are handled for any residue.** The grid code drops every non-backbone atom of the chosen residues, so the real side chain of an Arg or Gln does not block the path. The built-in surface check (`SurfaceTools.py`) uses the residue's actual atoms.
* **Result keys carry the residue name** as `(res_no, chain, res)`, and the output files write it.

## Proposed change

1. A new option, for example `-lys_sites FILE`, with one `chain:res` per line (author chain and residue number, as elsewhere). `jwalk` parses and validates the file.
2. After `read_PDB_file`, set `atom.res = "LYS"` for every atom whose `(chain, res_no)` is in the list, before the grid and CA marking. This does in memory what the user's converter does.
3. Report any listed site that is not found in the structure, as the crosslink-list code does for missing residues.
4. Add a test and a README entry for the new option.

Estimated size: about 30 to 50 lines, plus the test and documentation.

## The catch: `-surface` with FreeSASA

`check_solvent_accessibility_freesasa` (`SurfaceTools.py`) runs FreeSASA on the original file and matches its output to the marked residues using the residue name FreeSASA reports, which would still be ARG or GLN. With the in-memory rename the keys would not match, and those sites would silently count as buried.

The fix is to match on `(res_no, chain)` and take the residue name from the marked residues (about 5 lines). This is the part to test most carefully, with `-surface` on a structure where a listed site is not lysine.

## Open questions

* Should the output show the original residue name (ARG) or LYS? The rename approach shows LYS, the same as the user's converted PDB files do today.
* Should the override apply only with `-lys`, or with `-aa1`/`-aa2` too?
* List file format: `A:123` per line, or the pipe-delimited style of the crosslink list?

## Not done

Not implemented, and no code was run for this note. The line counts are estimates from reading the code. The work would go on a new branch, since `optimize` is already pushed.
