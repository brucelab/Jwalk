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

from array import array
import numpy as np
from multiprocessing import Pool, Process, freeze_support
import itertools
import sys
import math
import os

# order of voxels to search - by having diagonals last ensures shortest path is returned
COMB = [(1, 0, 0),
        (-1, 0, 0),
        (0, -1, 0),
        (0, 1, 0),
        (0, 0, -1),
        (0, 0, 1),
        (1, 0, 1),
        (-1, 0, 1),
        (0, 1, 1),
        (0, -1, 1),
        (1, -1, 0),
        (-1, -1, 0),
        (1, 1, 0),
        (-1, 1, 0),
        (1, 0, -1),
        (-1, 0, -1),
        (0, 1, -1),
        (0, -1, -1),
        (1, 1, 1),
        (1, -1, 1),
        (-1, 1, 1),
        (-1, -1, 1),
        (1, 1, -1),
        (1, -1, -1),
        (-1, 1, -1),
        (-1, -1, -1)]

def prepare_grid(dens_map, vox):

    '''

    Flattens the grid for the BFS. The grid is padded with a one voxel border that is
    marked as already visited, so the neighbours of a voxel sit at fixed offsets in the
    flat index and no bounds checks are needed.

    Returns (visited template, occupied voxels, y stride, x stride, steps) where steps is
    a list of (flat index offset, distance travelled) for each of the 26 neighbours.

    Arguments:

       *dens_map*
           grid with solvent accessible surface (masked array)
       *vox*
           number of angstoms per voxel

    '''

    nz, ny, nx = dens_map.fullMap.shape
    sy = ny + 2
    sx = nx + 2

    occupied = np.zeros((nz + 2, sy, sx), dtype=np.uint8)
    occupied[1:-1, 1:-1, 1:-1] = dens_map.fullMap > 0
    visited = np.ones((nz + 2, sy, sx), dtype=np.uint8)
    visited[1:-1, 1:-1, 1:-1] = 0

    # distance of diagonal steps
    diag1 = (math.sqrt((vox ** 2) * 2)) # 2d diagonal
    diag2 = (math.sqrt((vox ** 2) * 3)) # 3d diagonal
    cost = {1: vox, 2: diag1, 3: diag2}

    steps = []
    for (dx, dy, dz) in COMB:
        n_axes = (dx != 0) + (dy != 0) + (dz != 0)
        steps.append(((dz * sy + dy) * sx + dx, cost[n_axes]))

    return bytearray(visited.tobytes()), bytearray(occupied.tobytes()), sy, sx, steps

def breadth_first_search(start_voxels, grid, max_dist):

    '''

    Breadth First Search of grid from a set of starting voxels. For general info on
    algorithm see: https://en.wikipedia.org/wiki/Breadth-first_search

    Returns (visited, parent, root, distance) as flat arrays indexed like the grid made by
    prepare_grid. parent is the previous voxel on the path (-1 for a starting voxel), root
    is the starting voxel of the path and distance is the length of the path in angstroms.

    Arguments:

       *start_voxels*
           list of starting voxels [x, y, z]
       *grid*
           output of prepare_grid
       *max_dist*
           maximum distance BFS will search until

    '''

    visited_template, occupied, sy, sx, steps = grid

    visited = bytearray(visited_template) # voxels that have been searched
    n = len(visited)
    parent = array('i', [-1]) * n # voxel the path to this voxel came from
    root = array('i', [-1]) * n # starting voxel of the path to this voxel
    distance = array('d', [0.0]) * n # distance from starting voxel for each other voxel

    # place starting voxels into queue and initialise visited and root
    queue = []
    for (x, y, z) in start_voxels:
        p = ((z + 1) * sy + (y + 1)) * sx + (x + 1)
        visited[p] = 1
        root[p] = p
        queue.append(p)

    # grid is searched until queue is empty
    head = 0
    while head < len(queue):
        p = queue[head] # take first voxel in queue
        head += 1
        d = distance[p]
        if d <= max_dist:
            for (offset, step) in steps: # expand in all directions from voxel - in order of COMB.
                q = p + offset
                # check voxel hasn't already been searched (the border counts as searched)
                if not visited[q]:
                    visited[q] = 1
                    parent[q] = p
                    root[q] = root[p]
                    distance[q] = d + step
                    if not occupied[q]: # if the voxel is in empty space
                        queue.append(q) # add to queue for later searching

    return visited, parent, root, distance

def shortest_path(start_residue, end_residue, end_voxels, search, grid, aa1_CA, aa2_CA):

    '''

    Finds the shortest sasd to end_residue from the result of a BFS started at start_residue.

    Returns (length in angstroms, voxel path of sasd) or None if end_residue was not reached.

    Arguments:

       *start_residue*
           key of aa1_CA
       *end_residue*
           key of aa2_CA
       *end_voxels*
           list of ending voxels of end_residue
       *search*
           output of breadth_first_search
       *grid*
           output of prepare_grid
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha

    '''

    visited, parent, root, distance = search
    sy, sx = grid[2], grid[3]

    def unflatten(p):
        r, x = divmod(p, sx)
        z, y = divmod(r, sy)
        return [x - 1, y - 1, z - 1]

    (x_1, y_1, z_1) = aa1_CA[start_residue]
    (x_2, y_2, z_2) = aa2_CA[end_residue]

    shortest_distance = 9999
    shortest_voxel = None

    # cycling through possible end coords of end_residue to get shortest sasd
    for (x, y, z) in end_voxels:

        p = ((z + 1) * sy + (y + 1)) * sx + (x + 1)

        if visited[p]:

            (x_s, y_s, z_s) = unflatten(root[p])

            # add the distance between starting/ending residue CA voxel and start/end voxel in path
            d = distance[p]
            d += math.sqrt((x_1 - x_s) ** 2 + (y_1 - y_s) ** 2 + (z_1 - z_s) ** 2)
            d += math.sqrt((x - x_2) ** 2 + (y - y_2) ** 2 + (z - z_2) ** 2)

            # keep record of shortest distance - the last voxel wins a tie
            if d <= shortest_distance:
                shortest_distance = d
                shortest_voxel = p

    if shortest_voxel is None:
        return None

    # rebuild the path back to the starting voxel
    path = []
    p = shortest_voxel
    while p != -1:
        path.append(unflatten(p))
        p = parent[p]
    path.reverse()
    path.insert(0, aa1_CA[start_residue]) # add aa1 CA voxel to path
    path.append(aa2_CA[end_residue]) # add aa2 CA voxel to end of path

    return shortest_distance, path

def calculate_specific_SASD(single_crosslink, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                            max_dist, vox):

    '''

    Breadth First Search of grid. For general info on algorithm see:
    https://en.wikipedia.org/wiki/Breadth-first_search

    Returns dictionary containing solvent accessible surface distances between specific starting res
    and ending res.

    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *single_crosslink*
           start and end residue.
           start is key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel

    '''

    start_residue = single_crosslink[0]
    end_residue = single_crosslink[1]

    specific_xl = {}

    grid = prepare_grid(dens_map, vox)
    search = breadth_first_search(aa1_voxels[start_residue], grid, max_dist)

    # now we have a full set of paths into empty space starting from start_residue
    # all stored in search. Now need to extract paths to specific residue
    shortest = shortest_path(start_residue, end_residue, aa2_voxels[end_residue], search, grid,
                             aa1_CA, aa2_CA)

    if shortest:
        shortest_distance, path = shortest
        specific_xl[start_residue, end_residue, shortest_distance] = path # start lys, end lys, length of xl = path of xl

    return specific_xl


def calculate_SASDs(start_residue, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                    max_dist, vox):

    """

    Breadth First Search of grid. For general info on algorithm see:
    https://en.wikipedia.org/wiki/Breadth-first_search

    Returns dictionary containing solvent accessible surface distances between starting res
    and all possible ending res.

    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *start_residue*
           key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel


    """

    sasds = {}

    grid = prepare_grid(dens_map, vox)
    search = breadth_first_search(aa1_voxels[start_residue], grid, max_dist)

    # now we have a full set of paths into empty space starting from start_residue
    # all stored in search. Now need to extract paths to specific residues
    for end_residue in aa2_voxels:
        if start_residue != end_residue:

            shortest = shortest_path(start_residue, end_residue, aa2_voxels[end_residue], search,
                                     grid, aa1_CA, aa2_CA)

            # add shortest distance sasd to output dictionary

            if shortest:
                shortest_distance, path = shortest
                if start_residue[1] < end_residue[1]:  # this to order the dict so that chain goes alphabetically
                    sasds[start_residue, end_residue, shortest_distance] = path
                elif end_residue[1] < start_residue[1]:
                    sasds[end_residue, start_residue, shortest_distance] = path
                # if both on the same chain, then ordered to go numerically
                elif start_residue[0] < end_residue[0]:
                    sasds[start_residue, end_residue, shortest_distance] = path
                else:
                    sasds[end_residue, start_residue, shortest_distance] = path

    return sasds

def calculate_SASDs_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return calculate_SASDs(*a_b)
    
def calculate_specific_SASD_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return calculate_specific_SASD(*a_b)

def parallel_BFS(aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA, crosslink_pairs,
                 max_dist, vox, ncpus, xl_list):
    
    """
    
    Parallelised Breadth First Search of grid. 
    
    Returns dictionary containing all solvent accessible surface distances
    {start res, end res, length in angstroms : voxel path of sasd}
    
    Arguments:
    
       *start_residue*
           key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that 
           residue 
       *aa1_voxels*
           dictionary containing starting voxels {start_residue : starting voxels}
       *aa2_voxels*
           dictionary containing ending voxels {end_residue : ending voxels}
       *dens_map*
           grid with solvent accessible surface (masked array)
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha
       *crosslink_pairs*
           list of pairs of crosslinks (empty if not calculating specific crosslinks)
       *max_dist*
           maximum distance BFS will search until
       *vox*
           number of angstoms per voxel
       *ncpus*
           number of allocated cpus
           
    """ 
    
    freeze_support()
    final_XL = {}
    
    if xl_list:
        if ncpus > 1:
            pool = Pool(ncpus)
            xl_dictionaries = pool.map(calculate_specific_SASD_star,
                                       zip(crosslink_pairs,
                                       itertools.repeat(aa1_voxels),
                                       itertools.repeat(aa2_voxels),
                                       itertools.repeat(dens_map),
                                       itertools.repeat(aa1_CA),
                                       itertools.repeat(aa2_CA),
                                       itertools.repeat(max_dist),
                                       itertools.repeat(vox)))
                                       
            for c in xl_dictionaries:
                final_XL.update(c)

        else:
            # alternative call to allow single cpu running on Windows machines
            for single_crosslink in crosslink_pairs:
                xl_dictionaries = calculate_specific_SASD(single_crosslink, aa1_voxels,
                                                          aa2_voxels, dens_map, aa1_CA, aa2_CA,
                                                          max_dist, vox)
                final_XL.update(xl_dictionaries)

    else:
        if ncpus > 1:
            
            pool = Pool(ncpus)
            xl_dictionaries = pool.map(calculate_SASDs_star, 
                                       zip(aa1_voxels,
                                       itertools.repeat(aa1_voxels),
                                       itertools.repeat(aa2_voxels),
                                       itertools.repeat(dens_map),
                                       itertools.repeat(aa1_CA),
                                       itertools.repeat(aa2_CA),
                                       itertools.repeat(max_dist),
                                       itertools.repeat(vox)))
            
            for c in xl_dictionaries:
                final_XL.update(c)

        else:
            # alternative call to allow single cpu running on Windows machines
            for start_residue in aa1_voxels:
                xl_dictionaries = calculate_SASDs(start_residue, aa1_voxels, aa2_voxels,
                                                  dens_map, aa1_CA, aa2_CA, max_dist, vox)
                final_XL.update(xl_dictionaries)

    return final_XL

def calculate_distance(cords):
    ''' Calculates the distance of points in 3d, input e.g. [[x1,y1,z1],[x2,y2,z3]] '''
    return math.sqrt(((cords[0][0]-cords[1][0])**2)+((cords[0][1]-cords[1][1])**2)+((cords[0][2]-cords[1][2])**2))

def get_euclidean_distances(sasds, pdb, aa1, aa2):

    residues = {}
    euc_dists = {}
    with open (pdb) as inf:
        for line in inf:
            if line.startswith('ATOM') and (line[12:16].strip() == 'CA'):
                if line[21:22].strip() == "":
                    chain = " "
                else:
                    chain = line[21:22].strip()
                residues[line[22:26].strip(),chain] = [float(line[30:38].strip()),
                float(line[38:46].strip()),
                float(line[46:54].strip())]
    
    for k,v in residues.items():
        for k1,v1 in residues.items():
            if k1 != k:
            
                euc_dists[int(k[0]),k[1], int(k1[0]),k1[1]] = calculate_distance([v,v1])
    
    sasds_and_eucs = {}
    
    for s in sasds:
        if (s[0][0],s[0][1],s[1][0],s[1][1]) in euc_dists:
            sasds_and_eucs[s[0],s[1],s[2],euc_dists[(s[0][0],s[0][1],s[1][0],s[1][1])]] = sasds[s]
            
    return sasds_and_eucs

