#===============================================================================
#     This file is part of Jwalk.
#     
#     Jwalk - A tool to calculate the solvent accessible surface distance (SASD) 
#     between crosslinked residues.
#     
#     Copyright 2016 Jwalk Inventor and Birkbeck College University of London.
#                          The Jwalk Inventor is: Josh Bullock
# 
# 
#     Jwalk is available under Public Licence.
#     This software is made available under GPL V3
#
#     Please cite your use of Jwalk in published work:
#     
#     J.Bullock, J. Schwab, K. Thalassinos, M. Topf (2016)
#     The importance of non-accessible crosslinks and solvent accessible surface distance
#     in modelling proteins with restraints from crosslinking mass spectrometry. 
#     Molecular and Cellular Proteomics (15) pp.2491–2500
#
#===============================================================================

INSTALLATION:
Open a terminal session and a move to the directory Jwalk is in.

Create a virtual environment with the requirements (Python 3.8 or newer), from the
Jwalk directory:

  python3 -m venv src/venv
  src/venv/bin/pip install -r requirements.txt

then run src/jwalk-launcher in place of "jwalk" below. See docs/20261027_launcher.md.

RUNNING JWALK:

Jwalk currently only runs on Unix or Linux based operating systems (sorry !)

type "jwalk" plus any of the below flags:

  -h, --help        show this help message and exit
  -lys              calculate lysine crosslinks (default)
  -xl_list XL_LIST  calculate crosslinks from input list (see Examples)
  -i I              specify input structure, PDB or mmCIF: -i <inputfile.pdb> or -i <inputfile.cif>
                    (default runs on every .pdb, .cif and .mmcif file in the directory)
  -aa1 AA1          specify starting amino-acid via three letter code eg. HIS
  -aa2 AA2          specify ending amino-acid via three letter code eg. TYR (give both)
  -max_dist N       maximum crosslink distance in Angstroms, a whole number (default 60)
  -vox N            voxel size of the grid in Angstroms, a whole number of at least 1 (default 1)
  -ncpus N          number of cpus to use (default: all)
  -includebugs      reproduce the output of the original Jwalk (commit 72fac51): its path
                    extraction, which overestimates a SASD when two residues share a surface
                    voxel, and its surface spheres chosen by the first letter of the atom name
                    (so with -ligand, P and metal atoms get no sphere). For comparing with old
                    results; see docs/20261030-benchmark.md
  -surface          use higher accuracy method to calculate solvent accessibility - slower
                    needs the freesasa command on your PATH (FreeSASA 2.0 or 2.1,
                    http://freesasa.github.io/); with -ligand, freesasa also sees the ligands
  -ligand           treat non-water HETATM records (ligands, ions, cofactors) as obstacles
                    in the solvent accessible surface (default: HETATM ignored). With PDB input
                    this includes modified residues (MSE, ...); from mmCIF the modified residues
                    of a chain are always included, with or without -ligand
                    see docs/20261027-ligand-mode.md for what is included and how atoms
                    are sized

A value that isn't a whole number for -vox, -max_dist or -ncpus is reported as a usage error.

if no flags are added, Jwalk will calculate the SASDs between all lysine residues on
every .pdb, .cif and .mmcif file in the directory.

INPUT FORMATS:

PDB and mmCIF files are read alike (Biopython; chains and residue numbers are the author ones,
auth_asym_id and auth_seq_id, so mmCIF gives the same identifiers as the equivalent PDB). mmCIF
has no limit on the number of chains or the length of a chain id, and the modified residues of a
chain (MSE, PSU, OMG, ...), which are HETATM records, are always part of the structure because
the file says which residues belong to a chain. Ligands and ions still need -ligand. A PDB file
cannot tell a modified residue from a ligand, so there HETATM records are ligands. Results are
written under the file's base name (9syg.cif -> Jwalk_results/9syg_crosslink_list.txt).
See docs/20261030-mmcif-input.md.

Examples/test_ligand.pdb has a phosphate on the path between two lysines: run it with
and without -ligand to see the effect. tests/test_ligand.py checks ligand mode:

  src/venv/bin/python -m unittest discover -s tests -v

CROSSLINK LIST INPUT:

Crosslinks should be listed in a pipe-delimited .txt file as so:

res1|chain1|res2|chain2|

eg. 34|A|56|B|

If your pdb file has no chain, leave the chain fields empty

eg. 34||56||

(a lowercase x is looked up as a chain named x and finds nothing). mmCIF input always has
chains (auth_asym_id), so this applies to PDB files only.

See the example file for a sample list.

OUTPUT:

Two output files are generated in ./Jwalk_results/ 

<pdb>_crosslink_list.txt

<pdb>_crosslinks.pdb

In order to visualize the crosslinks in Chimera select Actions>Atoms/Bonds>show
Better representation in chimera in progress ...

Please try the test in JWalk/Examples/

Happy crosslinking !