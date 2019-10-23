"""
Created on 29/07/2013

@author: victor
"""
from __future__ import print_function


def symm_groups_validation( symm_groups):
    """
    Checks that symmetry groups are well defined (each n-tuple has a correspondent symmetric n-tuple)
    """
    try:
        for sg in symm_groups:
            for symm_pair in sg:
                if len(symm_pair) != 2:
                    raise Exception
    except Exception:
        raise ValueError('Symmetry groups are not well defined')



def symm_permutations(groups):
    """
    Generator that produces the possible permutations of the symmetry groups.
    Suposing we have this symmetries definition:
        [
            [ [1,2] ], 
            [ [3,4],[5,6] ]
        ]
    Meaning that coordinates 1 and 2 are equivalent, and that 3 and 4, as well as 5 and 6 are equivalen and must be changed 
    at the same time (they form a group). 
    The possible permutations are:
        [(1,2)] [(3,4) (5,6)]
        [(1,2)] [(4,3) (6,5)]
        [(2,1)] [(3,4) (5,6)]
        [(2,1)] [(4,3) (6,5)]
        
    :param groups: The symmetry groups list.
    
    :return: yields one permutation at a time.
    """
    if len(groups) > 0:
        head = groups[0]
        
        for tail_permutation in symm_permutations(groups[1:]):
            yield [head] + tail_permutation
        
        #swap all elements of the head group
        swapped_head = []
        for pair in head:
            swapped_head.append([pair[1],pair[0]])
        for tail_permutation in symm_permutations(groups[1:]):
            yield [swapped_head] + tail_permutation
    else:
        yield []
        
def swap_atoms(coordset_reference, atom_i, atom_j):
    #print("PRIMA", coordset_reference)
    coordset_reference[[atom_i,atom_j]] =  coordset_reference[[atom_j,atom_i]]
    #print("DOPO", coordset_reference)
    
def min_rmsd_of_rmsds_list(rmsds_list):
    return (rmsds_list.T).min(1)
    
