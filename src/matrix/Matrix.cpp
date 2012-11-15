#include "Statistics.h"
#include <python2.7/Python.h>
#include <python2.7/structmember.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

#define calc_vector_pos(i,j,matrix) (i*(matrix->row_length_minus_one) - i - ((( i-1)*i) >> 1 ) + j - 1)

typedef struct {
    PyObject_HEAD
    int row_length;
    int row_length_minus_one;
    int data_size;
    float* data;
    StatisticsCalculator* statisticsCalculator;
    bool statistics_already_calculated;
    PyObject* zero;
} CondensedMatrix;

/*
 * Object destructor. Only has to free memory used for data storage and the statistics calculator.
 */
static void condensedMatrix_dealloc(CondensedMatrix* self){

	// Special for this object
	delete [] self->data;

	// Python ops
	Py_XDECREF(self->zero);
    self->ob_type->tp_free((PyObject*)self);
}

/*
 *
 */
static PyObject* condensedMatrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	CondensedMatrix* self;
	self = (CondensedMatrix*) type->tp_alloc(type, 0);
    if (self != NULL) {
    	self->row_length = 0;
    	self->data_size = 0;
    	self->data = NULL;
    	self->zero =  Py_BuildValue("d", 0); // To be returned always that a 0 is needed
    	self->statisticsCalculator = NULL;
    	self->statistics_already_calculated = false;
    }
    return (PyObject *) self;
}

static int condensedMatrix_init(CondensedMatrix *self, PyObject *args, PyObject *kwds){
	PyObject* input = NULL;
	PyArrayObject* rmsd_numpy_array = NULL;
	bool numpy = false;
	if (! PyArg_ParseTuple(args, "O!",&PyArray_Type,&rmsd_numpy_array)){
		//cout<<"[CondensedMatrix] Getting matrix data from sequence."<<endl;
		if (! PyArg_ParseTuple(args, "O",&input)){
			PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters.");
			return -1;
		}

		rmsd_numpy_array = (PyArrayObject *) PyArray_ContiguousFromObject(input, PyArray_DOUBLE, 1, 1);
	}
	else{
		//cout<<"[CondensedMatrix] Getting matrix data from numpy array."<<endl;
		numpy = true;
	}

    if (rmsd_numpy_array == NULL){
    	PyErr_SetString(PyExc_RuntimeError, "Impossible to create intermediary data.\n"
    										"Check that the parameter is a sequence and there's memory available.");
    	return -1;
    }


    self->data_size = rmsd_numpy_array->dimensions[0];
    self->row_length = (int) (1 + sqrt(1+8*self->data_size))/2;
    self->row_length_minus_one = self->row_length - 1;
    self->data = new float[self->data_size];

    double * inner_data = (double*)(rmsd_numpy_array->data);

	for(int i = 0;i < self->data_size; ++i){
		self->data[i] = (float) (inner_data[i]);
	}

	// Let's alloc the statistics object
	self->statisticsCalculator =  new StatisticsCalculator(self->data,self->data_size);

	if(!numpy){
    	Py_DECREF(rmsd_numpy_array);
    }

	return 0;
}

static PyMemberDef condensedMatrix_members[] = {
    {"row_length", T_INT, offsetof(CondensedMatrix, row_length), READONLY,	PyDoc_STR("description")},
    {"data_size", T_INT, offsetof(CondensedMatrix, data_size), READONLY, PyDoc_STR("description")},
    {NULL}  /* Sentinel */
};

static PyObject* condensedMatrix_get_number_of_rows(CondensedMatrix* self, PyObject *args){
	return Py_BuildValue("i", self->row_length);
}

static PyObject* condensedMatrix_get_data(CondensedMatrix* self, PyObject *args){
	npy_intp dims[1] = {self->data_size};
	return  PyArray_SimpleNewFromData(1,dims,NPY_FLOAT,self->data);
}

static void condensedMatrix_calculate_statistics(CondensedMatrix* self, PyObject *args){
	self->statistics_already_calculated = false;
}

static PyObject* condensedMatrix_get_mean(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->mean);
}

static PyObject* condensedMatrix_get_variance(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->variance);
}

static PyObject* condensedMatrix_get_skewness(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->skewness);
}

static PyObject* condensedMatrix_get_kurtosis(CondensedMatrix* self, PyObject *args){
	if(self->statistics_already_calculated == false){
		self->statisticsCalculator->calculateStatistics();
		self->statistics_already_calculated = true;
	}
	return Py_BuildValue("d", self->statisticsCalculator->kurtosis);
}


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
	int inode,jnode,pos;
	double value;
	for (int i =0; i< len_nodes-1;++i){
		inode = nodes[i];
		for (int j = i+1; j< len_nodes;++j){
			jnode =nodes[j];
			pos = calc_vector_pos(inode,jnode,self);
			value =self->data[pos];
			if(value <= cutoff){
				neighbors[i] += 1;
				neighbors[j] += 1;
			}
		}
	}
	//Get index with maximum value
	int max =  *(max_element(neighbors.begin(),neighbors.end()));
	int index = 0;
	for (int i = 0; i< len_nodes-1;++i){
		if (neighbors[i]==max){
			index = i;
			break;
		}
	}

	PyObject* tuple = Py_BuildValue("(ii)", nodes[index],neighbors[index]);//PyTuple_Pack(2,nodes[index],neighbors[index]);
	delete [] nodes;
	return tuple;
}

/*
 def d(i,condensed_matrix):
    """
    Degree of a vertex:
    d_i = sum_{j=1}^n w_ij
    """
    d_val = 0
    n = condensed_matrix.row_length
    for j in range(n):
        d_val += condensed_matrix[i,j]
    return d_val
 */
double degree(int element, CondensedMatrix* self){
	double d_val = 0.;
	int pos;
	for (int i = 0;i<self->row_length;++i){
		if(element<i){
			pos = calc_vector_pos(element,i,self);
			d_val += (double) self->data[pos];
		}
		else{
			if (element > i){
				pos = calc_vector_pos(i,element,self);
				d_val += (double) self->data[pos];
			}
		}
	}
	return d_val;
}

static PyObject*  condensedMatrix_calculate_rw_laplacian(CondensedMatrix* self, PyObject *args){
	// First calculate D**-1, inverse diagonal matrix of degrees
	// (inverse of a diagonal matrix is the diagonal matrix in which Dinv(i,i) = 1/D(i,i)
	// W is indeed the rmsd matrix.
	// We shall calculate invD*W
	//	(d1       )     (w11 w21 w31 w41)    (d1*w11 d1*w21 d1*w31 d1*w41)
	//	(	d2    )  x  (w12 w22 w32 w42)  = (d2*w12 d2*w22 d2*w32 d2*w42)
	//	(	  d3  )     (w13 w23 w33 w43)    (d3*w13 d3*w23 d3*w33 d3*w43)
	//  (       d4)     (w14 w24 w34 w44)    (d4*w14 d4*w24 d4*w34 d4*w44)
	// Then to calculate I-Dinv*W only means to substract 1-Dinv_ij*W_ij for each i,j

	double invD;
	double * L = new double[self->row_length * self->row_length];
	double value;
	int pos;
	for (int i = 0; i<self->row_length; ++i){
		invD = 1. / degree(i,self);
		for (int j = 0; j<self->row_length; ++j){
			if(i==j){
				L[j + self->row_length*i] = 1.;
			}
			else{
				if( i<j){
					pos = calc_vector_pos(i,j,self);
				}
				else{
					pos = calc_vector_pos(j,i,self);
				}
				value =self->data[pos];
				L[j + self->row_length*i] = -1*value*invD;
			}
		}
	}
	//----------------------------
	// Is a transposition needed?
	//----------------------------
	npy_intp dims[2] = {self->row_length,self->row_length};
	return PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,L);
}
//dgeev
//dspevx
//dpteqr
static PyObject*  condensedMatrix_calculate_affinity_matrix(CondensedMatrix* self, PyObject *args){
	double delta;

	if (!PyArg_ParseTuple(args, "d",&delta)){
		PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters.");
		return NULL;
	}

	double value;
	int pos;
	double* affinity_matrix = new double[self->row_length * self->row_length];
	for (int i = 0; i<self->row_length; ++i){
			for (int j = 0; j<self->row_length; ++j){
				if(i==j){
					affinity_matrix[j + self->row_length*i] = 1.;
				}
				else{
					if( i<j){
						pos = calc_vector_pos(i,j,self);
					}
					else{
						pos = calc_vector_pos(j,i,self);
					}
					value =self->data[pos];
					affinity_matrix[j + self->row_length*i] = exp(-(value*value)/(2*delta*delta));
				}
			}
		}


	npy_intp dims[2] = {self->row_length,self->row_length};
	return PyArray_SimpleNewFromData(2,dims,NPY_DOUBLE,affinity_matrix);
}

static PyMethodDef condensedMatrix_methods[] = {
	// Basic field access
	{"get_number_of_rows", (PyCFunction) condensedMatrix_get_number_of_rows, METH_NOARGS,PyDoc_STR("description")},
	{"get_data", (PyCFunction) condensedMatrix_get_data, METH_NOARGS,PyDoc_STR("description")},

	// Statistics
	{"recalculateStatistics", (PyCFunction) condensedMatrix_calculate_statistics, METH_NOARGS,PyDoc_STR("description")},
	{"calculateMean", 		(PyCFunction) condensedMatrix_get_mean, METH_NOARGS,PyDoc_STR("description")},
	{"calculateVariance", 	(PyCFunction) condensedMatrix_get_variance, METH_NOARGS,PyDoc_STR("description")},
	{"calculateSkewness", 	(PyCFunction) condensedMatrix_get_skewness, METH_NOARGS,PyDoc_STR("description")},
	{"calculateKurtosis", 	(PyCFunction) condensedMatrix_get_kurtosis, METH_NOARGS,PyDoc_STR("description")},

	// Matrix as graph
	{"get_neighbors_for_node", (PyCFunction)condensedMatrix_get_neighbors_for_node, METH_VARARGS,PyDoc_STR("description")},
	{"choose_node_with_higher_cardinality", (PyCFunction)condensedMatrix_choose_node_with_higher_cardinality, METH_VARARGS,PyDoc_STR("description")},
	{"calculate_rw_laplacian", (PyCFunction)condensedMatrix_calculate_rw_laplacian, METH_NOARGS,PyDoc_STR("description")},
	{"calculate_affinity_matrix", (PyCFunction)condensedMatrix_calculate_affinity_matrix, METH_VARARGS,PyDoc_STR("description")},
	{NULL}  /* Sentinel */
};

PyObject* condensedMatrix_subscript(CondensedMatrix *self, PyObject *key){
	int pos, i,j;
	i = PyInt_AS_LONG(PyTuple_GET_ITEM(key,0));
	j = PyInt_AS_LONG(PyTuple_GET_ITEM(key,1));

	if (i < j){
		pos = calc_vector_pos(i,j,self);
		return PyFloat_FromDouble((double)self->data[pos]);
	}
	else{
		if (i==j){
			Py_INCREF(self->zero);
			return self->zero;
		}
		else{
			pos = calc_vector_pos(j,i,self);
			return PyFloat_FromDouble((double)self->data[pos]);
		}
	}
}

int condensedMatrix_ass_subscript(CondensedMatrix *self, PyObject *key, PyObject *v){
	int pos, i,j;
	i = PyInt_AS_LONG(PyTuple_GET_ITEM(key,0));
	j = PyInt_AS_LONG(PyTuple_GET_ITEM(key,1));

	if (i < j){
		pos = calc_vector_pos(i,j,self);
		self->data[pos] = PyFloat_AsDouble(v);
	}
	else{
		if (i!=j){
			pos = calc_vector_pos(j,i,self);
			self->data[pos] = PyFloat_AsDouble(v);
		}
	}
	return 0;
}

Py_ssize_t condensedMatrix_length(CondensedMatrix *self){
	return self->row_length;
}

static PyMappingMethods pdb_as_mapping = {
    (lenfunc)     	condensedMatrix_length,				/*mp_length*/
    (binaryfunc)	condensedMatrix_subscript,			/*mp_subscript*/
    (objobjargproc)	condensedMatrix_ass_subscript,		/*mp_ass_subscript*/
};

static PyTypeObject CondensedMatrixType = {
    PyObject_HEAD_INIT(NULL)
    0,                         						/*ob_size*/
    "condensedMatrix.CondensedMatrixType",      	/*tp_name*/
    sizeof(CondensedMatrix), 	/*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)condensedMatrix_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    &pdb_as_mapping,           /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT| Py_TPFLAGS_BASETYPE , /*tp_flags*/
    "Condensed matrix as in pdist",           /* tp_doc */
	0,		               					  /* tp_traverse */
	0,		               					  /* tp_clear */
	0,		               					  /* tp_richcompare */
	0,		              					  /* tp_weaklistoffset */
	0,		               	   /* tp_iter */
	0,		               	   /* tp_iternext */
	condensedMatrix_methods,   /* tp_methods */
	condensedMatrix_members,   /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)condensedMatrix_init, /* tp_init */
	0,                         		/* tp_alloc */
	condensedMatrix_new,        		/* tp_new */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif


PyMODINIT_FUNC initcondensedMatrix(void){
    PyObject* module;

    if (PyType_Ready(&CondensedMatrixType) < 0)
        return;

    module = Py_InitModule3("condensedMatrix", NULL,"Fast Access Condensed Matrix");
    if (module == NULL)
          return;

    Py_INCREF(&CondensedMatrixType);
    PyModule_AddObject(module, "CondensedMatrix", (PyObject*) &CondensedMatrixType);

    import_array();
}
