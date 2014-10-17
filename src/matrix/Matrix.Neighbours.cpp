/*
 def get_neighbors_for_node(condensed_matrix,node,nodes_left,cutoff):
    """
    As the name of the function says, it will pick all the neighbor elements for a given
    element of the dataset. One element is neighbor of another if their distance falls within
    a cutoff. nodes_lef will be the remaining nodes of the dataset, so we'll find the neighbors
    there.
    """
    neighbours = []
    # Scan the column
    for i in range(len(nodes_left)):
        if node > nodes_left[i]:
            #print nodes_left[i], node
            # hoping there's lazy evaluation...
            if condensed_matrix[nodes_left[i],node]<= cutoff:#access_element(condensed_matrix, nodes_left[i], node-1, row_len) <= cutoff:
                neighbours.append(nodes_left[i])
    # Scan the row
    for i in range(len(nodes_left)):
        if node < nodes_left[i]:
            #print  node, nodes_left[i]
            if condensed_matrix[node, nodes_left[i]]<= cutoff: #access_element(condensed_matrix, node, nodes_left[i]-1, row_len) <= cutoff:
                neighbours.append(nodes_left[i])
    return neighbours

 */

static PyObject* condensedMatrix_get_neighbors_for_node(CondensedMatrix* self, PyObject *args){
	// Parse all arguments
	PyObject *nodes_left_list, *node_obj,*cutoff_obj;

	if (!PyArg_ParseTuple(args, "OOO",&node_obj,&nodes_left_list,&cutoff_obj)){
		PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters.");
		return NULL;
	}

	// Convert to C types
	int node = (int) PyInt_AsLong((PyObject *)node_obj);
	double cutoff = PyFloat_AsDouble((PyObject *)cutoff_obj);
	int len_nodes_left = PyList_Size((PyObject *)nodes_left_list);
	int* nodes_left = new int[len_nodes_left];
	for(int i = 0; i < len_nodes_left; ++i){
		nodes_left[i] = PyInt_AS_LONG(PyList_GET_ITEM((PyObject*) nodes_left_list,i));
	}

	// Do the job
	vector<int> neighbours;
	int pos,j;
	for(int i = 0; i < len_nodes_left; ++i){
		j = nodes_left[i];
		if(node<j){
			pos = calc_vector_pos(node,j,self);
			if(self->data[pos]<=cutoff){
				neighbours.push_back(j);
			}
		}
		if(node>j){
			pos = calc_vector_pos(j,node,self);
			if(self->data[pos]<=cutoff){
				neighbours.push_back(j);
			}
		}
	}

	int neigh_len = neighbours.size();
	npy_intp dims[1] = {neigh_len};
	int* neighbours_data = new int[neigh_len];
	copy(&(neighbours[0]),&(neighbours[0]) + neigh_len,neighbours_data);
	delete [] nodes_left;
	return PyArray_SimpleNewFromData(1,dims,NPY_INT,neighbours_data);
}
/*
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
*/
static PyObject*  condensedMatrix_choose_node_with_higher_cardinality(CondensedMatrix* self, PyObject *args){
	// Parse all arguments
	PyObject *nodes_list,*cutoff_obj;

	if (!PyArg_ParseTuple(args, "OO",&nodes_list,&cutoff_obj)){
		PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters.");
		return NULL;
	}

	// Convert to C types
	double cutoff = PyFloat_AsDouble((PyObject *)cutoff_obj);
	int len_nodes = PyList_Size((PyObject *)nodes_list);
	int* nodes = new int[len_nodes];
	for(int i = 0; i < len_nodes; ++i){
		nodes[i] = PyInt_AS_LONG(PyList_GET_ITEM((PyObject*) nodes_list,i));
	}
	
	//Do the job
	vector<int> neighbors(len_nodes,0);
	int inode, jnode, pos;
	double value;
	for (int i = 0; i< len_nodes-1; ++i){
		inode = nodes[i];
		for (int j = i+1; j < len_nodes; ++j){
			jnode = nodes[j];
			pos = calc_vector_pos(inode, jnode, self);
			value = self->data[pos];
			if(value <= cutoff){
				neighbors[i] += 1;
				neighbors[j] += 1;
			}
		}
	}
	
	//Get maximum value
	int max =  *(max_element(neighbors.begin(),neighbors.end()));
	// Get the index that corresponds to that value
	// cout<<"DBG: Looking for "<<max<<endl;
	int index = 0;
	for (int i = 0; i < len_nodes; ++i){
		if (neighbors[i]==max){
			index = i;
			break;
		}
	}
//	cout<<"DBG: (node, nneig) ";
//    for (int i = 0; i < len_nodes; ++i){
//        cout<<"("<<nodes[i]<<", "<<neighbors[i]<<") ";
//    }
//    cout<<endl;
//    cout<<"DBG: nodes, neig "<<nodes[index]<<" "<<neighbors[index]<<endl;
	PyObject* tuple = Py_BuildValue("(ii)", nodes[index],neighbors[index]);//PyTuple_Pack(2,nodes[index],neighbors[index]);
	delete [] nodes;
	return tuple;
}

/*
 # __eps_neighborhood will not return the central element!!
def __eps_neighborhood(self,node,eps):
	"""
	Populates the neighbor list for a given node (without the node itself).
	"""
	neighbours = []
	for i in range(self.number_of_elements):
		if node != i and self.condensed_matrix[node,i] <= eps:
			neighbours.append(i)
	return neighbours

 */
// Simpler version of the neighbouring algorithm
static PyObject* condensedMatrix_get_neighbors_of_node_for_radius(CondensedMatrix* self, PyObject *args){
	// Parse all arguments
	double radius = 0.0;
	int node = -1;
	long int number_of_nodes = self->row_length;

	if (!PyArg_ParseTuple(args, "id",&node, &radius)){
		PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters.");
		return NULL;
	}

	// Do the job
	vector<int> neighbours;
	int pos;
	for(int i = 0; i < number_of_nodes; ++i){
		if(i < node){
			pos = calc_vector_pos(i, node, self);
			if(self->data[pos] <= radius){
				neighbours.push_back(i);
			}
		}

		if(i > node){
			pos = calc_vector_pos(node, i, self);
			if(self->data[pos] <= radius){
				neighbours.push_back(i);
			}
		}

		//if i == node, do nothing
	}

	int neigh_len = neighbours.size();
	npy_intp dims[1] = {neigh_len};
	int* neighbours_data = new int[neigh_len];
	copy(&(neighbours[0]), &(neighbours[0]) + neigh_len, neighbours_data);
	return PyArray_SimpleNewFromData(1,dims,NPY_INT,neighbours_data);
}
