'''
Created on 29/07/2013

@author: victor
'''

def symm_group_permutator(symm_groups_left, used, permutations):
   
    if len(symm_groups_left) > 0:
        # 2 options: not to use it
        symm_group_permutator(symm_groups_left[1:],used, permutations)
        # Or using it
        used_now = list(used)
        used_now.append(symm_groups_left[0])
        symm_group_permutator(symm_groups_left[1:], used_now, permutations)
    else:
        permutations.append(list(used))
        
def swap_atoms(coordset_reference, atom_i, atom_j):
    coordset_reference[[atom_i,atom_j]] =  coordset_reference[[atom_j,atom_i]]