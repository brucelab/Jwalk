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
#     Molecular and Cellular Proteomics (15) pp.2491-2500
#
#===============================================================================

import sys
from math import cos, sin
import subprocess
from Jwalk import PDBTools
import os
import re
import string

def expand_points(atom,sphere):
    """ Exapnd the unit sphere around specific x,y,z point and returns surface points"""
    points = {}
    atom.point_list = []
    CH2 = 1.68
    radius = {"N":1.43,
              "O":1.30,
              "C":1.68,
              "S":1.67,
              }
    
    r = radius[atom.atom_name[0]] + CH2
    (x,y,z) = (atom.x, atom.y, atom.z)
    
    for s in sphere:
        x1 = x + s[0]*r
        y1 = y + s[1]*r
        z1 = z + s[2]*r
        points[x1,y1,z1] = 0

    return points

def create_unit_sphere():
    """ Generates a unit sphere with 30 points on the surface """
    
    unit_sphere = []
    
    unit_sphere.append( [0.0, 1.0, 0.0] )
    unit_sphere.append( [0.0, -1.0, 0.0] )
    
    nstep = 5
    PI = 3.1415926536
    theta = PI/nstep
    arc = theta
        
    for istep in range(nstep):
        istep = istep + 1    # to change range from 0--9 to 1--10
        y1 = cos(istep*theta)
        r2 = sin(istep*theta)
        ndot2= 2*PI*r2/arc   # the circumference at that radius / proportion of pi
        if ndot2 == 0.0:
            continue
        theta2 = 2*PI/ndot2
        for idot in range(int(ndot2)):
            idot = idot + 1  # to change range from 0-- to 1--
            x2 = r2*cos(idot*theta2)
            z2 = r2*sin(idot*theta2)
            unit_sphere.append( [x2, y1, z2] )
    
    return unit_sphere 
	
def check_solvent_accessibility(prot,aa1_CA,xl_list = False):
    
    '''
    
    Checks solvent accessibility of residues in aa1_CA
    Returns aa1_CA of solvent accessible residues. 
    
    Arguments
        *prot*
            Tempy structure instance
        *aa1_CA*
            residues of interest
            
    '''
    # dictionary for text output
    string_dict = {"LYS":"lysines",
                "CYS":"cysteines",
                "ASP":"acidic residues",
                "GLU":"acidic residues", 
                "VAL":"valines",
                "ILE":"isoleucines",
                "LEU":"leucines",
                "ARG":"arginines",
                "PRO":"prolines",
                "GLY":"glycines",
                "ALA":"alanines",
                "TRP":"tryptophans",
                "PHE":"phenylalanines",
                "SER":"serines",
                "GLN":"glutamines",
                "HIS":"histidines",
                "MET":"methionines",
                "THR":"threonines",
                "ASN":"asparagines",
                "TYR":"tyrosines"
                }
    
    radius = {"N":1.43,
              "O":1.30,
              "C":1.68,
              "S":1.67,
              }
    sd_res = False
    # create sphere of 30 points to expand around atoms
    sphere = create_unit_sphere()
    
    SA_res = {}
    CH2 = 1.68
    
    # this is not very efficient ... freesasa implementation to come
    for atom in prot.atomList:
        atom.clash = 1
        if (atom.res_no,atom.chain,atom.res) in aa1_CA:
            sd_res = atom.res
            # generate 30 points in unit spehere around atom
            points = expand_points(atom,sphere)
            
            for p in points:
                (x,y,z) = (p[0],p[1],p[2])
                # for every other atom check if points intersects with it
                for atom2 in prot.atomList:
                    if atom2.res != atom.res or (atom2.res == atom.res and atom2.atom_name != atom.atom_name):
                        r = radius[atom2.atom_name[0]] + CH2
                        # need to transpose x,y,z
                        (tx,ty,tz) = (x-atom2.x,y-atom2.y,z-atom2.z)
                        # if the point lies within the sphere of that atom then it clashes
                        if tx**2 + ty**2 + tz**2 <= r**2:
                            points[p] = 1
                            break
                # if any point on the sphere doesn't intersect with another then the atom is solvent accessible
                if points[p] == 0:
                    atom.clash = 0
                    break
            
            # if atom doesn't clash then residue information is kept in SA_res
            if atom.clash == 0:
                SA_res[atom.res_no,atom.chain,atom.res] = aa1_CA[atom.res_no,atom.chain,atom.res]
    
    # inform user on buried resiudes
    if xl_list:
        pass
    elif sd_res == "LYS":
        print("%d %s and 1 N-terminus of which %d are on the surface" % (len(aa1_CA)-1,string_dict[sd_res], len(SA_res)))                
    else:            
        print("%d %s of which %d are on the surface" % (len(aa1_CA),string_dict[sd_res], len(SA_res)))
            
    return SA_res   
    
def update_crosslink_pairs(crosslink_pairs, aa1_CA, aa2_CA, remove_aa1, remove_aa2):
    
    '''Removes buried residues from crosslink_pairs'''
    
    index_to_delete = set()
    
    # sets, so checking each residue pair against the residues to remove is quick
    remove_aa1_set = set((x[0], x[1]) for x in remove_aa1)
    remove_aa2_set = set((x[0], x[1]) for x in remove_aa2)
    
    for i in range(len(crosslink_pairs)): # for each residue pair, check both are solvent accessible
        
        x1, x2 = crosslink_pairs[i]
        
        if x1 not in aa1_CA:
            index_to_delete.add(i)
        elif x2 not in aa2_CA:
            index_to_delete.add(i)
                
        if (x1[0],x1[1]) in remove_aa1_set:
            index_to_delete.add(i)
        elif (x2[0],x2[1]) in remove_aa2_set:        
            index_to_delete.add(i)
    
    no_sasd_possible = []
    crosslink_pairs_final = []
    for i in range(len(crosslink_pairs)):
        if i not in index_to_delete:
            crosslink_pairs_final.append(crosslink_pairs[i])
        else:
            no_sasd_possible.append(crosslink_pairs[i])
    
    if len(no_sasd_possible) > 0:
        print("the following crosslinks cannot be calculated:")
        for s in no_sasd_possible:
            print("%s-%s-%s - %s-%s-%s" % (s[0][2],s[0][0],s[0][1],s[1][2],s[1][0],s[1][1]))
                    
    return crosslink_pairs_final
    
def rsa_words(residues):

    """

    Returns {(residue name, chain + residue number): (res_no, chain, res)} for the residues of
    interest, the form in which a freesasa .rsa made from an mmCIF file names a residue.

    Arguments:

       *residues*
           {(res_no, chain, res) : ...}

    """

    return dict(((res, chain.strip() + str(num)), (num, chain, res)) for (num, chain, res) in residues)

def parse_rsa_line(line, words=None):

    """

    Reads one RES line of a freesasa .rsa file.

    Returns (residue number, chain, residue name, relative all-atom SASA) or None if the line is
    not a residue line or has no relative SASA (freesasa writes N/A for hetero residues).

    The columns are separated by whitespace rather than fixed, because freesasa 2.1 moved the
    chain to make room for longer chain labels. With a chain the line is
    RES LYS A 21 ABS REL ... (14 fields), without one RES LYS 21 ABS REL ... (13 fields).
    An insertion code on the residue number is dropped, as it is everywhere else in Jwalk.

    With mmCIF input freesasa runs the chain label and the residue number together (A022 is
    residue 22 of chain A0). Such a line is matched against words, the output of rsa_words for
    the residues of interest, first as written and then without a trailing insertion code.

    Arguments:

       *line*
           one line of the .rsa file
       *words*
           output of rsa_words, for a .rsa made from an mmCIF file

    """

    if not line.startswith("RES"):
        return None

    col = line.split()

    if words and len(col) == 13:
        word = col[2]
        key = (col[1], word) if (col[1], word) in words else (col[1], word.rstrip(string.ascii_letters))
        if key in words:
            num, chain, name = words[key]
            try:
                return num, chain, name, float(col[4])
            except ValueError:
                return None

    if len(col) == 14:
        name, chain, num, rel = col[1], col[2], col[3], col[5]
    elif len(col) == 13:
        name, chain, num, rel = col[1], " ", col[2], col[4]
    else:
        return None

    number = re.match(r"-?\d+", num)

    try:
        return int(number.group()), chain, name, float(rel)
    except (AttributeError, ValueError):
        return None

def write_pdb_without_water(prot, filename):

    """

    Copies a pdb file leaving out the water molecules, so that freesasa --hetatm sees the same
    atoms as Jwalk does with -ligand (read_PDB_file keeps hetero residues but not water; it
    treats HOH and WAT as water, as Biopython does).

    Arguments:

       *prot*
           .pdb file to copy
       *filename*
           file to write

    """

    with open(prot) as inf, open(filename, "w") as outf:
        for line in inf:
            if line.startswith("HETATM") and line[17:20].strip() in ("HOH", "WAT"):
                continue
            outf.write(line)

def check_solvent_accessibility_freesasa(prot,aa1_CA,freesasa_source = "freesasa", xl_list = False, hetatm = False):

    """

    Runs freesasa on the pdb file and returns the residues of aa1_CA that it finds on the surface
    (relative all-atom SASA above 7%).

    Arguments:

       *prot*
           .pdb or .cif file
       *aa1_CA*
           residues of interest {(res_no, chain, res) : CA voxel}
       *freesasa_source*
           freesasa executable
       *xl_list*
           True if the residues come from a crosslink list
       *hetatm*
           include hetero residues (ligands, ions) in the freesasa calculation, as -ligand does
           in the Jwalk grid. Water is left out either way.

    """

    # dictionary for text output
    string_dict = {"LYS":"lysines",
                "CYS":"cysteines",
                "ASP":"aspartates",
                "GLU":"glutamates", 
                "VAL":"valines",
                "ILE":"isoleucines",
                "LEU":"leucines",
                "ARG":"arginines",
                "PRO":"prolines",
                "GLY":"glycines",
                "ALA":"alanines",
                "TRP":"tryptophans",
                "PHE":"phenylalanines",
                "SER":"serines",
                "GLN":"glutamines",
                "HIS":"histidines",
                "MET":"methionines",
                "THR":"threonines",
                "ASN":"asparagines",
                "TYR":"tyrosines"
                }
    
    rsa_file = os.path.splitext(prot)[0] + ".rsa"

    # an .rsa left by an earlier run must not be mistaken for this one's
    if os.path.exists(rsa_file):
        os.remove(rsa_file)

    cmd = [freesasa_source,"--format=rsa","--output=%s" % rsa_file,"--radii=naccess"]

    mmcif = PDBTools.is_mmcif(prot)

    if mmcif:
        # freesasa reads mmCIF with --cif. It is given a copy holding exactly the atoms of Jwalk's
        # grid (the chains with their modified residues, ligands only with -ligand, no water) and
        # --hetatm so that it counts the HETATM rows that are left.
        input_pdb = os.path.splitext(prot)[0] + "_jwalk.cif"
        PDBTools.write_mmcif_for_surface(prot, input_pdb, hetatm)
        cmd += ["--cif", "--hetatm", input_pdb]
    elif hetatm:
        # freesasa --hetatm would also count water, which -ligand leaves out of the grid
        input_pdb = os.path.splitext(prot)[0] + "_nowater.pdb"
        write_pdb_without_water(prot, input_pdb)
        cmd += ["--hetatm", input_pdb]
    else:
        input_pdb = None
        cmd.append(prot)

    try:
        ext_freesasa = subprocess.Popen(cmd)
        ext_freesasa.communicate()
    except OSError:
        print("Could not run '%s'. Please check you have Freesasa installed and on your PATH. Go to 'http://freesasa.github.io/' " % freesasa_source)
        exit(1)
    finally:
        if input_pdb and os.path.exists(input_pdb):
            os.remove(input_pdb)

    if not os.path.exists(rsa_file):
        print("No .rsa file created. Please check you have Freesasa installed. Go to 'http://freesasa.github.io/' ")
        exit(1)

    # for an mmCIF file freesasa writes the chain label and residue number as one word (A022 for
    # residue 22 of chain A0), so the residues of interest are looked up by that word
    words = rsa_words(aa1_CA) if mmcif else {}

    SA_res = {}
    with open(rsa_file) as rsa:
        for line in rsa:
            parsed = parse_rsa_line(line, words)
            if parsed:
                NUM, CHAIN, RES, REL = parsed
                if REL > 7.0:
                    SA_res[(NUM,CHAIN,RES)] = True

    SA_res_dict = {}
    
    for s in aa1_CA:
        if (s[0],s[1],s[2]) in SA_res:
            SA_res_dict[(s[0],s[1],s[2])] = aa1_CA[(s[0],s[1],s[2])]
            sd_res = s[2]
        else:
            print("Residue %d-%s-%s is buried" % (s[0],s[1],s[2]))
            sd_res = s[2]
    
    # inform user on buried resiudes
    if xl_list:
        pass
    elif sd_res == "LYS":
        print("%d %s and 1 N-terminus of which %d are on the surface" % (len(aa1_CA)-1,string_dict[sd_res], len(SA_res_dict)))                
    else:            
        print("%d %s of which %d are on the surface" % (len(aa1_CA),string_dict[sd_res], len(SA_res_dict)))
            
    return SA_res_dict
    