"""
Connectity patterns used in Basal Ganglia modeling literature.

@author     Lucas Koelman
"""

from __future__ import division
from enum import unique
import numpy as np
from   scipy.linalg import circulant

from bgcellmodels.common.stdutil import IntEnumDescriptor


def adjacency_to_list(adj_mat, src_dim=0, threshold=0):
    """
    Convert adjacency matrix to list connection list.

    @param      adj_mat : np.array
                Two-dimensional adjacency matrix

    @param      src_dim : int
                Dimension of adjacency matrix that corresponds to source
                population.

    @param      threshold : float
                Threshold value for connections in adjacency matrix.

    @return     conn_list : list(list<int, int>)
                List of pairs [i, j] where i is the cell index in source population
                and j the cell index in target population
    """
    if src_dim == 1:
        adj_mat = adj_mat.T # transpose
    elif src_dim != 0:
        raise ValueError("Source dimensions must be either '0' or '1'.")
    if adj_mat.ndim != 2:
        raise ValueError("Adjacency matrix must be two-dimensional array.")
    
    src_target_ids = np.where(adj_mat > threshold)
    conn_list = zip(*src_target_ids)
    return conn_list


@unique
class ConnectivityPattern(IntEnumDescriptor):
    """
    Named connectivity patterns between STN and GPe populations
    """
    # 1-to-1 intra-channel
    # (Kumaravelu 2016: THA-RCX, RCX-DSTR, RCX_ISTR)
    Kumaravelu2016_OneToOne_Intra = 1

    # 2-to-1 (1 intra-channel, 1 adjacent)
    # (Kumaravelu 2016: GPE-STN, STN-GPE, RCX-STN, GPE-GPI)
    Kumaravelu2016_TwoToOne_IntraAdj = 2

    # 2-to-1 on some neurons, 0-to-1 on others
    # (Kumaravelu 2016: STN-GPE, STN-GPI)
    Kumaravelu2016_TwoToSome_StnGpi = 3  # see comment below
    # Fig, 1: 'some GPe/GPi neurons receive excitatory input from two STN neurons, while others do not'.
    # MATLAB code: project to intra-channel and right adjacent on 20%/50% of GPe/GPi neurons. 
    # This means 50% of neurons receive afferent from intra-channel and left adjacent.
    Kumaravelu2016_TwoToSome_StnGpe = 4
    Kumaravelu2016_GpeGpe = 9

    # Article So (2012)
    So2012_GpeStn = 5
    So2012_GpeGpi = 6
    So2012_GpeGpe = 7
    So2012_FromSTN = 8       # fig. 1: 'Each STN cell projects to two neighbouring GPe and GPi cells'

    # Rubin & Terman (2002) and (2004) articles
    RubinTerman_RandomSparse_StnGpe = 10       # fig. 3a in Rubin & Terman (2002)
    RubinTerman_RandomSparse_GpeStn = 11       # fig. 3a in Rubin & Terman (2002)
    RubinTerman_RandomSparse_GpeGpe = 12       # fig. 3a in Rubin & Terman (2002)

    RubinTerman_StructuredSparse_StnGpe = 13   # fig. 5a in Rubin & Terman (2002)
    RubinTerman_StructuredSparse_GpeStn = 14   # fig. 5a in Rubin & Terman (2002)
    RubinTerman_StructuredSparse_GpeGpe = 15   # fig. 5a in Rubin & Terman (2002)
    
    RubinTerman_StructuredTight_StnGpe = 16    # fig. 7a in Rubin & Terman (2002)
    RubinTerman_StructuredTight_GpeStn = 17    # fig. 7a in Rubin & Terman (2002)
    RubinTerman_StructuredTight_GpeGpe = 18    # fig. 7a in Rubin & Terman (2002)

    RubinTerman2004_StnGpe = 19    # approximation of Rubin & Terman (2004)
    RubinTerman2004_GpeStn = 20    # approximation of Rubin & Terman (2004)
    RubinTerman2004_GpeGpe = 21    # approximation of Rubin & Terman (2004)


def make_divergent_pattern(src, tar_rel, pop_size):
    """
    Make fixed divergent connection pattern with wrap-around (periodic
    boundaries).
    
    @param      src : list(int)
                Indices of source cells, e.g. range(num_cell) for all cells
    
    @param      tar_rel : list(int)
                Relative indices of target cells, e.g. [0, -1]

    @return     conn_list : list(list<int, int>)
                List of pairs [i, j] where i is the cell index in source population
                and j the cell index in target population
    """
    return [[i, (i+j)%pop_size] for i in src for j in tar_rel]


def make_convergent_pattern(tar, src_rel, pop_size):
    """
    Make fixed convergent connection pattern with wrap-around (periodic
    boundaries).

    @param      tar : list(int)
                Indices of target cells, e.g. range(num_cell) for all cells
    
    @param      src_rel : list(int)
                Relative indices of source cells, e.g. [0, -1]

    @return     conn_list : list(list<int, int>)
                List of pairs [i, j] where i is the cell index in source population
                and j the cell index in target population
    """
    return [[(i+j)%pop_size, i] for i in tar for j in src_rel]


def make_divergent_pattern_rt2004(num_cell, num_target, sep, offset, rng):
    """
    Make divergent connectivity pattern based on a generalization of 
    the connectivity pattern in the code of Rubin & Terman (2004).

    @param      num_cell : int
                Number of cells per population (number of 'channels')

    @param      num_target : int OR list([int, int])
                Number of cells targeted by one pre-synaptic cell: fixed number
                or range.

    @param      sep : list([int, int])
                Range for separation between target cells (in number of cells),
                i.e. the value for j_b-j_a if j is a target cell index and a
                and b two neighboring target cells.

    @param      offset : list([int, int])
                Range for the offset from the pre-synaptic cell index of the 
                middle target cell (or last in bottom half if num_target is even).

    @param      rng : numpy.random
                Numpy random generator
    """
    conn_list = []

    # Number of target cells for each pre-synaptic cell
    if not isinstance(num_target, int):
        num_target = rng.randint(num_target[0], num_target[1]+1, num_cell)
    else:
        num_target = [num_target for i in range(num_cell)]
    
    for i in range(num_cell):
        # separation between target cels is between [8, 2]
        n_target = num_target[i]
        separation = rng.randint(sep[0], sep[1]+1, n_target-1)
        # offset of middle target is between [-3, 3]
        j_offset, = rng.randint(offset[0], offset[1]+1, 1)
        j_start = i + j_offset - sum(separation[0:int(len(separation)/2)])
        if j_start < 0:
            j_start = 0
        overshoot = j_start + sum(separation) - (num_cell-1)
        if overshoot > 0:
            j_start -= overshoot
        # make target cell indices
        for jj in range(n_target):
            j_target = j_start + sum(separation[0:jj])
            conn_list.append([i, j_target])

    return conn_list


def make_connection_list(pattern, num_cell, min_pre=None, rng=None):
    """
    Make connection lists based on named connectivity pattern.

    @param      pattern : ConnectivityPattern
                Enum instance representing connectivity pattern

    @param      num_cell : int
                Number of cells per population (number of 'channels')

    @param      rng : numpy.random
                Random number generator (seeded)

    @return     conn_list : list(list<int, int>)
                List of pairs [i, j] where i is the cell index in source population
                and j the cell index in target population
    """
    if rng is None:
        rng = np.random

    pattern = ConnectivityPattern.get(pattern) # accept strings

    #===========================================================================
    # Kumaravelu (2016)
    if pattern == ConnectivityPattern.Kumaravelu2016_OneToOne_Intra:
        conn_list = [[i,i] for i in range(num_cell)]
    
    elif pattern == ConnectivityPattern.Kumaravelu2016_TwoToOne_IntraAdj:
        tar = range(num_cell)
        src_rel = [0,1] # relative channel numbers
        conn_list = make_convergent_pattern(tar, src_rel, num_cell)
    
    elif pattern == ConnectivityPattern.Kumaravelu2016_TwoToSome_StnGpi:
        tar = rng.choice(range(num_cell), int(0.5*num_cell), 
                               replace=False) # pick 50% of target layer randomly
        src_rel = [0,-1]
        conn_list = make_convergent_pattern(tar, src_rel, num_cell)
    
    elif pattern == ConnectivityPattern.Kumaravelu2016_TwoToSome_StnGpe:
        # See comment above
        tar = rng.choice(range(num_cell), int(0.2*num_cell), replace=False) # pick 20% of target layer randomly
        src_rel = [0,-1]
        conn_list = make_convergent_pattern(tar, src_rel, num_cell)
    
    elif pattern == ConnectivityPattern.Kumaravelu2016_GpeGpe:
        # to left adjacent channel and second right neighbour
        src = range(num_cell)
        tar_rel = [-1,2]
        conn_list = make_divergent_pattern(src, tar_rel, num_cell)
    
    #===========================================================================
    # So (2012)
    elif pattern == ConnectivityPattern.So2012_FromSTN:
        src = range(num_cell)
        tar_rel = [0,1]
        conn_list = make_divergent_pattern(src, tar_rel, num_cell)
    
    elif pattern == ConnectivityPattern.So2012_GpeGpe:
        src = range(num_cell)
        tar_rel = [-1,2]
        conn_list = make_divergent_pattern(src, tar_rel, num_cell)
    
    elif pattern == ConnectivityPattern.So2012_GpeGpi: # same as GpeGpe
        src = range(num_cell)
        tar_rel = [-1,2]
        conn_list = make_divergent_pattern(src, tar_rel, num_cell)
    
    elif pattern == ConnectivityPattern.So2012_GpeStn:
        src = range(num_cell)
        tar_rel = [-1,0]
        conn_list = make_divergent_pattern(src, tar_rel, num_cell)
    
    #===========================================================================
    # Rubin & Terman (2002)
    elif pattern == ConnectivityPattern.RubinTerman_RandomSparse_StnGpe:
        # 1 STN cell -> 1 random GPe cell
        js = rng.permutation(num_cell)
        conn_list = [[i,j] for i,j in enumerate(js)]

    elif pattern == ConnectivityPattern.RubinTerman_RandomSparse_GpeStn:
        # 1 GPe cell -> 3 random STN cells
        conn_list = [[i,j] for i in range(num_cell) for j in rng.choice(num_cell, 
                                                            3, replace=False)]

    elif pattern == ConnectivityPattern.RubinTerman_RandomSparse_GpeGpe:
        # 1 GPe cell -> 6 neighbouring GPe cells
        neighbours = [-3,-2,-1,1,2,3]
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]
    
    elif pattern == ConnectivityPattern.RubinTerman_StructuredSparse_StnGpe:
        # 1 STN cell -> 1 specific GPe cell (topographic/intra-channel)
        conn_list = [[i,i] for i in range(num_cell)]

    elif pattern == ConnectivityPattern.RubinTerman_StructuredSparse_GpeStn:
        # 1 GPe cell -> 2 specific STN cells (2nd neighbouring channel)
        neighbours = [-2,2]
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]

    elif pattern == ConnectivityPattern.RubinTerman_StructuredSparse_GpeGpe:
        # 1 GPe cell -> 2 specific GPe cell (2 nearest neighbors)
        neighbours = [-1,1]
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]
    
    elif pattern == ConnectivityPattern.RubinTerman_StructuredTight_StnGpe:
        # 1 STN cell -> 3 specific GPe cells (topographic, same and neighbouring channels)
        neighbours = [-1,0,1]
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]

    elif pattern == ConnectivityPattern.RubinTerman_StructuredTight_GpeStn:
        # 1 GPe cell -> 5 specific STN cells (topographic, same and neighbouring channels)
        neighbours = range(-2,3)
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]

    elif pattern == ConnectivityPattern.RubinTerman_StructuredTight_GpeGpe:
        # 1 GPe cell -> 6 specific GPe cells (neighbours)
        neighbours = [-3,-2,-1,1,2,3]
        conn_list = [[i,(i+j)%num_cell] for i in range(num_cell) for j in neighbours]
    
    #===========================================================================
    # Rubin & Terman (2004)
    elif pattern == ConnectivityPattern.RubinTerman2004_StnGpe:
        # 1 STN cell -> 3 GPe cells separated by 2-4 cells, topographic, with wrap-around
        # conn_list = [
        #     [0,2],[0,7],[0,10], [1,3],[1,6],[1,11], [2,1],[2,5],[2,13], 
        #     [3,0],[3,4],[3,12], [4,2],[4,6],[4,14], [5,3],[5,7],[5,15], 
        #     [6,1],[6,5],[6,9], [7,0],[7,4],[7,8], [8,2],[8,10],[8,15],
        #     [9,3],[9,11],[9,14], [10,5],[10,9],[10,13], [11,4],[11,8],[11,12],
        #     [12,6],[12,10],[12,14], [13,7],[13,11],[13,15], [14,1],[14,9],[14,13],
        #     [15,0],[15,8],[15,12]]
        conn_list = make_divergent_pattern_rt2004(num_cell, 3, [2, 8], [-3, 3], rng)

    elif pattern == ConnectivityPattern.RubinTerman2004_GpeStn:
        # 1 GPe cell -> 2 STN cells separated by 2-4 cells, topographic, no wrap-around
        # conn_list = [
        #     [0,1],[0,5], [1,0],[1,4], [2,3],[2,6], [3,2],[3,7], [4,0],[4,5],
        #     [5,1],[5,4], [6,3],[6,7], [7,2],[7,6], [8,9],[8,13], [9,8],[9,12],
        #     [10,11],[10,14], [11,10],[11,15], [12,8],[12,13], [13,9],[13,12],
        #     [14,11],[14,15], [15,10],[15,14]]
        conn_list = make_divergent_pattern_rt2004(num_cell, 2, [3, 5], [-4, 1], rng)

    elif pattern == ConnectivityPattern.RubinTerman2004_GpeGpe:
        # 1 GPe cell -> 1-3 GPe cells separated by 1-4 cells, topographic, no wrap-around
        # conn_list = [
        #     [0,1],[0,3], [1,0],[1,5], [2,0],[2,3],[2,6], [3,2],[3,7], 
        #     [4,1],[4,5], [5,4], [6,4],[6,7], [7,2],[7,6], [8,9],[8,11],
        #     [9,8],[9,13], [10,8],[10,11],[10,14], [11,10],[11,15],
        #     [12,9],[12,13], [13,12], [14,12],[14,15], [15,10],[15,14]]
        conn_list = make_divergent_pattern_rt2004(num_cell, [1, 3], [3, 5], [-5, -1], rng)

    else:
        ValueError('Unknown connection pattern <{0}>'.format(pattern))

    # Duplicate projections if minimum number of inputs is requested for each
    # post-synaptic cell.
    if min_pre is not None:
        for j in range(num_cell):
            inputs = [proj for proj in conn_list if proj[1] == j]
            num_inputs = len(inputs)
            num_extra = min_pre - num_inputs
            while num_inputs > 0 and num_extra > 0:
                # keep duplicating until we have enough inputs
                conn_list.append(inputs[num_extra % num_inputs])
                num_extra -= 1
    
    return conn_list

###############################################################################
# Watts-Strogatz Small-World network
###############################################################################

# Implementation Copyright Francis Song, see http://www.nervouscomputer.com/hfs/super-simple-watts-strogatz/
 
def _distance_matrix(L):
    Dmax = L//2
 
    D  = range(Dmax+1)
    D += D[-2+(L%2):0:-1]
 
    return circulant(D)/Dmax
 
def _pd(d, p0, beta):
    return beta*p0 + (d <= p0)*(1-beta)
 

def watts_strogatz(L, p0, beta, directed=False, rngseed=1):
    """
    Watts-Strogatz model of a small-world network
 
    This generates the full adjacency matrix, which is not a good way to store
    things if the network is sparse.

    License
    -------

    Copyright Francis Song, see http://www.nervouscomputer.com/hfs/super-simple-watts-strogatz/
 

    Parameters
    ----------
    L        : int
               Number of nodes.
 
    p0       : float
               Edge density. If K is the average degree then p0 = K/(L-1).
               For directed networks "degree" means out- or in-degree.
 
    beta     : float
               "Rewiring probability."
 
    directed : bool
               Whether the network is directed or undirected.
 
    rngseed  : int
               Seed for the random number generator.
 

    Returns
    -------
    A        : (L, L) array
               Adjacency matrix of a WS (potentially) small-world network.
    """
    rng = np.random.RandomState(rngseed)
 
    d = _distance_matrix(L)
    p = _pd(d, p0, beta)
 
    if directed:
        A = 1*(rng.random_sample(p.shape) < p)
        np.fill_diagonal(A, 0)
    else:
        upper = np.triu_indices(L, 1)
 
        A          = np.zeros_like(p, dtype=int)
        A[upper]   = 1*(rng.rand(len(upper[0])) < p[upper])
        A.T[upper] = A[upper]
 
    return A


def small_world_list(L, p0, beta, directed=False, rngseed=1):
    """
    Watts-Strogatz small-world network as connection list.

    @see    watts_strogatz() for meaning of parameters
    """
    adj_mat = watts_strogatz(L, p0, beta, directed, rngseed)
    return adjacency_to_list(adj_mat)