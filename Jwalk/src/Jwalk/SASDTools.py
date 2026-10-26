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
from collections import OrderedDict
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

    Prepares the grid for the BFS.

    Returns (occupied voxels, costs, vox) where occupied is an array of the grid that is 1 for
    voxels inside the solvent accessible surface and costs gives the distance travelled for each
    of the 26 neighbours, in the order of COMB.

    Arguments:

       *dens_map*
           grid with solvent accessible surface (masked array)
       *vox*
           number of angstoms per voxel

    '''

    # distance of diagonal steps
    diag1 = (math.sqrt((vox ** 2) * 2)) # 2d diagonal
    diag2 = (math.sqrt((vox ** 2) * 3)) # 3d diagonal
    cost = {1: vox, 2: diag1, 3: diag2}

    costs = [cost[(dx != 0) + (dy != 0) + (dz != 0)] for (dx, dy, dz) in COMB]

    return (dens_map.fullMap > 0).astype(np.uint8), np.array(costs, dtype=np.float64), vox

def prepare_targets(voxel_lists):

    '''

    Collects the voxels a search is looking for, so the search can stop once it has found them all.

    Returns an array of the target voxels [x, y, z] without repeats.

    Arguments:

       *voxel_lists*
           iterable of lists of target voxels [x, y, z]

    '''

    voxels = [voxel for voxel_list in voxel_lists for voxel in voxel_list]

    return np.unique(np.array(voxels, dtype=np.int64).reshape(-1, 3), axis=0)

def breadth_first_search(start_voxels, grid, max_dist, targets=None):

    '''

    Breadth First Search of grid from a set of starting voxels. For general info on
    algorithm see: https://en.wikipedia.org/wiki/Breadth-first_search

    The search goes a level at a time. Candidate voxels are laid out in the order a queue would
    reach them (each voxel in the level, then each of its neighbours in the order of COMB) and
    the first time a voxel appears wins, so the paths and distances found are the same as for a
    voxel at a time search with a queue.

    A step along a path is never shorter than the number of voxels it moves along any axis, so
    the search can't get further than max_dist / vox voxels (plus the step that finds a voxel)
    from the starting voxels. Only that part of the grid, padded with a one voxel border that is
    marked as already searched, is held in memory, and the neighbours of a voxel sit at fixed
    offsets in the flat index so no bounds checks are needed.

    Returns (visited, parent, root, distance, window). visited, parent, root and distance are
    flat arrays over the window. parent is the previous voxel on the path (-1 for a starting
    voxel), root is the starting voxel of the path and distance is the length of the path in
    angstroms. window is (x, y, z of the first voxel, then x, y, z size) of the part of the grid
    held, and is used by shortest_path to find a voxel in the arrays.

    Arguments:

       *start_voxels*
           list of starting voxels [x, y, z]
       *grid*
           output of prepare_grid
       *max_dist*
           maximum distance BFS will search until
       *targets*
           output of prepare_targets. The search stops as soon as all of these voxels have been
           reached - the path to a voxel is fixed when it is first reached, so this doesn't
           change any result for them.

    '''

    occupied_grid, costs, vox = grid
    n_moves = len(costs)
    grid_size = np.array(occupied_grid.shape[::-1]) # x, y, z

    start = np.array(start_voxels, dtype=np.int64).reshape(-1, 3)

    # the part of the grid the search can reach
    if len(start) > 0:
        margin = int(math.ceil(max_dist / vox)) + 2
        first = np.maximum(start.min(axis=0) - margin, 0)
        last = np.minimum(start.max(axis=0) + margin + 1, grid_size)
    else:
        first = last = np.zeros(3, dtype=np.int64)
    x_0, y_0, z_0 = (int(i) for i in first)
    x_size, y_size, z_size = (int(i) for i in last - first)
    sx = x_size + 2
    sy = y_size + 2
    n = (z_size + 2) * sy * sx

    inside = (slice(1, -1),) * 3

    visited = np.ones((z_size + 2, sy, sx), dtype=np.uint8) # voxels that have been searched
    visited[inside] = 0
    visited = visited.reshape(-1)
    occupied = np.zeros((z_size + 2, sy, sx), dtype=np.uint8)
    occupied[inside] = occupied_grid[z_0:z_0 + z_size, y_0:y_0 + y_size, x_0:x_0 + x_size]
    occupied = occupied.reshape(-1)

    offsets = np.array([(dz * sy + dy) * sx + dx for (dx, dy, dz) in COMB], dtype=np.int64)

    index = np.int32 if n < 2 ** 31 - 1 else np.int64
    parent = np.full(n, -1, dtype=index) # voxel the path to this voxel came from
    root = np.full(n, -1, dtype=index) # starting voxel of the path to this voxel
    distance = np.zeros(n, dtype=np.float64) # distance from starting voxel for each other voxel

    unclaimed = np.iinfo(np.int64).max
    claim = np.full(n, unclaimed, dtype=np.int64) # scratch space to find the first path to a voxel

    window = (x_0, y_0, z_0, x_size, y_size, z_size)

    # place starting voxels into the first level and initialise visited and root
    level = ((start[:, 2] - z_0 + 1) * sy + (start[:, 1] - y_0 + 1)) * sx + (start[:, 0] - x_0 + 1)
    visited[level] = 1
    root[level] = level

    if targets is not None:
        # flag the target voxels that are in the part of the grid that is held
        found = targets - first
        found = found[((found >= 0) & (found < last - first)).all(axis=1)]
        flags = np.zeros(n, dtype=np.uint8)
        flags[((found[:, 2] + 1) * sy + (found[:, 1] + 1)) * sx + (found[:, 0] + 1)] = 1
        remaining = int(flags.sum()) - int(flags[np.unique(level)].sum())
        if remaining <= 0:
            return visited, parent, root, distance, window

    # grid is searched until there are no more voxels to search from
    while len(level):

        # voxels further than max_dist are not searched from
        level = level[distance[level] <= max_dist]
        if len(level) == 0:
            break

        # expand in all directions from every voxel in the level - in order of COMB.
        candidates = (level[:, None] + offsets).ravel()

        # check voxel hasn't already been searched (the border counts as searched)
        unsearched = np.flatnonzero(visited[candidates] == 0)
        if len(unsearched) == 0:
            break

        # a voxel reached more than once belongs to the first path that reaches it
        found = candidates[unsearched]
        position = np.arange(len(found))
        np.minimum.at(claim, found, position) # lowest position that reaches each voxel
        first_reached = np.flatnonzero(claim[found] == position)
        claim[found] = unclaimed
        new = found[first_reached]
        source = unsearched[first_reached]
        previous = level[source // n_moves]

        visited[new] = 1
        parent[new] = previous
        root[new] = root[previous]
        distance[new] = distance[previous] + costs[source % n_moves]

        if targets is not None: # stop once every target voxel has been reached
            remaining -= int(flags[new].sum())
            if remaining <= 0:
                break

        level = new[occupied[new] == 0] # if the voxel is in empty space it is searched from

    return visited, parent, root, distance, window

def shortest_path(start_residue, end_residue, end_voxels, search, aa1_CA, aa2_CA):

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
       *aa1_CA*
           dictionary containing voxel of C-alpha
       *aa2_CA*
           dictionary containing voxel of C-alpha

    '''

    visited, parent, root, distance, window = search
    x_0, y_0, z_0, x_size, y_size, z_size = window
    sx = x_size + 2
    sy = y_size + 2

    def unflatten(p):
        r, x = divmod(p, sx)
        z, y = divmod(r, sy)
        return [x - 1 + x_0, y - 1 + y_0, z - 1 + z_0]

    (x_1, y_1, z_1) = aa1_CA[start_residue]
    (x_2, y_2, z_2) = aa2_CA[end_residue]

    shortest_distance = 9999
    shortest_voxel = None

    # cycling through possible end coords of end_residue to get shortest sasd
    for (x, y, z) in end_voxels:

        # a voxel outside the part of the grid that was held can't have been reached
        if not (x_0 <= x < x_0 + x_size and y_0 <= y < y_0 + y_size and z_0 <= z < z_0 + z_size):
            continue

        p = ((z - z_0 + 1) * sy + (y - y_0 + 1)) * sx + (x - x_0 + 1)

        if visited[p]:

            (x_s, y_s, z_s) = unflatten(int(root[p]))

            # add the distance between starting/ending residue CA voxel and start/end voxel in path
            d = float(distance[p])
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
        p = int(parent[p])
    path.reverse()
    path.insert(0, aa1_CA[start_residue]) # add aa1 CA voxel to path
    path.append(aa2_CA[end_residue]) # add aa2 CA voxel to end of path

    return shortest_distance, path

def calculate_specific_SASDs(start_residue, end_residues, aa1_voxels, aa2_voxels, dens_map,
                             aa1_CA, aa2_CA, max_dist, vox, grid=None):

    '''

    Breadth First Search of grid. For general info on algorithm see:
    https://en.wikipedia.org/wiki/Breadth-first_search

    One search from start_residue is shared by all of its specific crosslinks.

    Returns a list with a dictionary for each end residue in end_residues, in the same order,
    containing the solvent accessible surface distance between start_residue and that end residue

    {start res, end res, length in angstroms : voxel path of sasd}

    Arguments:

       *start_residue*
           key of aa1_voxels. aa1_voxels[start_residue] = all the starting voxels for that
           residue
       *end_residues*
           list of keys of aa2_voxels, the residues start_residue is crosslinked to
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
       *grid*
           output of prepare_grid, made from dens_map and vox if not given

    '''

    if grid is None:
        grid = prepare_grid(dens_map, vox)

    targets = prepare_targets([aa2_voxels[end_residue] for end_residue in end_residues])
    search = breadth_first_search(aa1_voxels[start_residue], grid, max_dist, targets)

    # now we have the paths into empty space starting from start_residue all stored in search.
    # Now need to extract paths to specific residues
    specific_xls = []
    for end_residue in end_residues:

        specific_xl = {}

        shortest = shortest_path(start_residue, end_residue, aa2_voxels[end_residue], search,
                                 aa1_CA, aa2_CA)

        if shortest:
            shortest_distance, path = shortest
            specific_xl[start_residue, end_residue, shortest_distance] = path # start lys, end lys, length of xl = path of xl

        specific_xls.append(specific_xl)

    return specific_xls

def calculate_specific_SASD(single_crosslink, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                            max_dist, vox, grid=None):

    '''

    Returns dictionary containing solvent accessible surface distances between specific starting res
    and ending res. See calculate_specific_SASDs.

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
       *grid*
           output of prepare_grid, made from dens_map and vox if not given

    '''

    return calculate_specific_SASDs(single_crosslink[0], [single_crosslink[1]], aa1_voxels,
                                    aa2_voxels, dens_map, aa1_CA, aa2_CA, max_dist, vox, grid)[0]


def calculate_SASDs(start_residue, aa1_voxels, aa2_voxels, dens_map, aa1_CA, aa2_CA,
                    max_dist, vox, grid=None, targets=None):

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
       *grid*
           output of prepare_grid, made from dens_map and vox if not given
       *targets*
           output of prepare_targets for all of aa2_voxels, made if not given


    """

    sasds = {}

    if grid is None:
        grid = prepare_grid(dens_map, vox)
    if targets is None:
        targets = prepare_targets(aa2_voxels.values())
    search = breadth_first_search(aa1_voxels[start_residue], grid, max_dist, targets)

    # now we have the paths into empty space starting from start_residue, all stored in search.
    # Now need to extract paths to specific residues
    for end_residue in aa2_voxels:
        if start_residue != end_residue:

            shortest = shortest_path(start_residue, end_residue, aa2_voxels[end_residue], search,
                                     aa1_CA, aa2_CA)

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

# data shared by every task, set once per worker process so it is not pickled with each task
shared = {}

def init_worker(grid, targets, aa1_voxels, aa2_voxels, aa1_CA, aa2_CA, max_dist):
    """Store the data shared by every task in the worker process."""
    shared.update(grid=grid, targets=targets, aa1_voxels=aa1_voxels, aa2_voxels=aa2_voxels,
                  aa1_CA=aa1_CA, aa2_CA=aa2_CA, max_dist=max_dist)

def calculate_SASDs_task(start_residue):
    """Run calculate_SASDs for one start residue using the shared data."""
    return calculate_SASDs(start_residue, shared['aa1_voxels'], shared['aa2_voxels'], None,
                           shared['aa1_CA'], shared['aa2_CA'], shared['max_dist'], None,
                           grid=shared['grid'], targets=shared['targets'])

def calculate_specific_SASDs_task(crosslinks):
    """Run calculate_specific_SASDs for a start residue and its end residues using the shared data."""
    start_residue, end_residues = crosslinks
    return calculate_specific_SASDs(start_residue, end_residues, shared['aa1_voxels'],
                                    shared['aa2_voxels'], None, shared['aa1_CA'], shared['aa2_CA'],
                                    shared['max_dist'], None, grid=shared['grid'])

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

    # the grid is flattened once and shared by every search
    grid = prepare_grid(dens_map, vox)

    if xl_list:
        # crosslinks that start at the same residue share a search
        task = calculate_specific_SASDs_task
        end_residues = OrderedDict()
        pair_numbers = OrderedDict()
        for i, (start_residue, end_residue) in enumerate(crosslink_pairs):
            end_residues.setdefault(start_residue, []).append(end_residue)
            pair_numbers.setdefault(start_residue, []).append(i)
        residues = list(end_residues.items())
        targets = None
    else:
        task = calculate_SASDs_task
        residues = list(aa1_voxels)
        targets = prepare_targets(aa2_voxels.values())

    initargs = (grid, targets, aa1_voxels, aa2_voxels, aa1_CA, aa2_CA, max_dist)

    if ncpus > 1:
        # results come back in the order of residues so the output order doesn't depend on timing
        with Pool(ncpus, init_worker, initargs) as pool:
            xl_dictionaries = list(pool.imap(task, residues, chunksize=1))

    else:
        # alternative call to allow single cpu running on Windows machines
        init_worker(*initargs)
        xl_dictionaries = [task(residue) for residue in residues]

    if xl_list:
        # put the crosslinks back in the order they were given
        by_pair = [None] * len(crosslink_pairs)
        for start_residue, c in zip(pair_numbers, xl_dictionaries):
            for i, specific_xl in zip(pair_numbers[start_residue], c):
                by_pair[i] = specific_xl
        xl_dictionaries = by_pair

    for c in xl_dictionaries:
        final_XL.update(c)

    return final_XL

def calculate_distance(cords):
    ''' Calculates the distance of points in 3d, input e.g. [[x1,y1,z1],[x2,y2,z3]] '''
    return math.sqrt(((cords[0][0]-cords[1][0])**2)+((cords[0][1]-cords[1][1])**2)+((cords[0][2]-cords[1][2])**2))

def get_euclidean_distances(sasds, pdb, aa1, aa2):

    residues = {}
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
    
    coords = {}
    for k,v in residues.items():
        coords[int(k[0]),k[1]] = v

    sasds_and_eucs = {}
    
    # only the distances between the ends of each sasd are needed
    for s in sasds:
        start = (s[0][0],s[0][1])
        end = (s[1][0],s[1][1])
        if start != end and start in coords and end in coords:
            sasds_and_eucs[s[0],s[1],s[2],calculate_distance([coords[start],coords[end]])] = sasds[s]
            
    return sasds_and_eucs

