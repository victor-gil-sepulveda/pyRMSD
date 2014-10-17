"""
Created on 16/11/2012

@author: victor
"""
import numpy

def choose_node_with_higher_cardinality(condensed_matrix,nodes,cutoff):
    """
    Returns the node in 'nodes' which has the bigger number of neighbours. One node
    is a neighbour of other if the distance between them is lower than the cutoff.
    The distances are stored in a condensed form distance matrix ('condensed_matrix') which
    represents a 'row_len' x 'row_len' symmetric square matrix.
    """
    neighbors = numpy.array([0]*len(nodes))
    len_nodes = len(nodes)
    nodes.sort()
    for i in range(len_nodes-1):
        inode = nodes[i]
        for j in range(i+1,len_nodes):
            #print inode, nodes[j],":",access_element(condensed_matrix, inode, nodes[j]-1, row_len)
            if condensed_matrix[inode,nodes[j]]<=cutoff: #access_element(condensed_matrix, inode, nodes[j]-1, row_len) <= cutoff:
                neighbors[i] += 1
                neighbors[j] += 1
        #print neighbors
    idx = neighbors.argmax()
    return (nodes[idx],neighbors[idx])

def get_neighbors_for_node(condensed_matrix,node,nodes_left,cutoff):
    """
    As the name of the function says, it will pick all the neighbor elements for a given
    element of the dataset. One element is neighbor of another if their distance falls within
    a cutoff. nodes_lef will be the remaining nodes of the dataset, so we'll find the neighbors
    there.
    """
    neighbours = []
    #print "columns"
    # Scan the column
    for i in range(len(nodes_left)):
        if node > nodes_left[i]: 
            #print nodes_left[i], node
            # hoping there's lazy evaluation...
            if condensed_matrix[nodes_left[i],node]<= cutoff:#access_element(condensed_matrix, nodes_left[i], node-1, row_len) <= cutoff:
                neighbours.append(nodes_left[i])
    #print "row"
    # Scan the row
    for i in range(len(nodes_left)):
        if node < nodes_left[i]:
            #print  node, nodes_left[i]
            if condensed_matrix[node, nodes_left[i]]<= cutoff: #access_element(condensed_matrix, node, nodes_left[i]-1, row_len) <= cutoff:
                neighbours.append(nodes_left[i])
    return neighbours