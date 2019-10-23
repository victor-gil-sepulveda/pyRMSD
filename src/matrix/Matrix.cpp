#include "Statistics.h"
#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
using namespace std;

#define calc_vector_pos(i,j,matrix) (i*(matrix->row_length_minus_one) - i - ((( i-1)*i) >> 1 ) + j - 1)

#define InputType int
#define INPUT_TYPE_NONE 0
#define INPUT_TYPE_LIST 1
#define INPUT_TYPE_NUMPY 2

#if PY_VERSION_HEX >= 0x03000000
#define PyInt_AsLong(x) PyLong_AsLong(x)
#define PyInt_AS_LONG PyLong_AS_LONG
#endif



typedef struct {
    PyObject_HEAD
    InputType input_type;
    long int row_length;
    long int data_size;
    // Data
    float* data;
    PyObject* numpy_array;
    // Statistics
    StatisticsCalculator* statisticsCalculator;
    bool statistics_already_calculated;
    // Precalculated stuff
    PyObject* zero;
    long int row_length_minus_one;
} CondensedMatrix;

/*
 * Object destructor. Only has to free memory used for data storage and the statistics calculator.
 */
static void condensedMatrix_dealloc(CondensedMatrix* self){

	// Special for this object
	delete [] self->data;
	delete self->statisticsCalculator;

	// Python ops
	Py_XDECREF(self->zero);
	if(self->numpy_array != NULL){
		Py_DECREF(self->numpy_array);
	}
    Py_TYPE(self)->tp_free((PyObject*)self);
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
    	self->numpy_array = NULL;
    	self->zero =  Py_BuildValue("d", 0.); // To be returned always that a 0 is needed
    	self->statisticsCalculator = NULL;
    	self->statistics_already_calculated = false;
    	self->input_type = INPUT_TYPE_NONE;
    }
    return (PyObject *) self;
}

static int condensedMatrix_init(CondensedMatrix *self, PyObject *args, PyObject *kwds){
	PyObject* input = NULL;

	if (!PyArg_ParseTuple(args, "O",&input)){
		PyErr_SetString(PyExc_ValueError, "[CondensedMatrix] Undefined problem parsing parameters.");
		return -1;
	}

	if (PyArray_Check(input)){
		//cout<<"[CondensedMatrix] Getting matrix data from numpy array."<<endl;
		self->input_type = INPUT_TYPE_NUMPY;
	}
	else{
		if (PyList_Check(input)){
			self->input_type = INPUT_TYPE_LIST;
			//cout<<"[CondensedMatrix] Getting matrix data from list."<<endl;
		}
		else{
			self->input_type = INPUT_TYPE_NONE;
		}

	}

	switch(self->input_type){
		case INPUT_TYPE_LIST:
		{
			self->data_size = PyList_Size(input);

			if (self->data_size < 1 || !PyFloat_Check(PyList_GetItem(input, 0))){
				PyErr_SetString(PyExc_ValueError, "[CondensedMatrix] Input list must be a 1D real vector with size > 0.");
				return -1;
			}
			self->row_length = (long int) (1 + sqrt(1+8*self->data_size))/2;
			self->row_length_minus_one = self->row_length - 1;
			self->data = new float[self->data_size];
			for(int i = 0;i < self->data_size; ++i){
				self->data[i] = (float) PyFloat_AS_DOUBLE(PyList_GetItem(input, i));
			}

			npy_intp dims[] = {(npy_intp)self->data_size};
			self->numpy_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, self->data);
		}
		break;

		case INPUT_TYPE_NUMPY:
		{
			//cout<<"[CondensedMatrix] Processing numpy array."<<endl;
			PyArrayObject* numpy_array = (PyArrayObject*) input;

			if (numpy_array->descr->type_num != NPY_DOUBLE || numpy_array->nd != 1)  {
				PyErr_SetString(PyExc_ValueError, "[CondensedMatrix] Input numpy array must be a 1D real vector.");
				return -1;
			}

			self->data_size = PyArray_DIMS(numpy_array)[0];
			self->row_length = (long int) (1 + sqrt(1+8*self->data_size))/2;
			self->row_length_minus_one = self->row_length - 1;
			double * inner_data = (double*)  PyArray_GETPTR1(numpy_array,0);
			self->data = new float[self->data_size];

			for(int i = 0; i < self->data_size; ++i){
				self->data[i] = (float) (inner_data[i]);
			}

			self->numpy_array = input;
			Py_INCREF(self->numpy_array);
		}
		break;

		case INPUT_TYPE_NONE:
			PyErr_SetString(PyExc_RuntimeError, "[CondensedMatrix] Input must be either a list or a numpy array.");
			return -1;
			break;

		default:
			PyErr_SetString(PyExc_RuntimeError, "[CondensedMatrix] Unexpected input behaviour.");
			return -1;
	}


	// Let's alloc the statistics object
	self->statisticsCalculator =  new StatisticsCalculator(self->data, self->data_size);
	return 0;
}

char row_length_text[] = "row_length";
char data_size_text[] = "data_size";
char bogus_description_text[] = "TODO";

static PyMemberDef condensedMatrix_members[] = {
    {row_length_text, T_INT, offsetof(CondensedMatrix, row_length), READONLY,	PyDoc_STR(bogus_description_text)},
    {data_size_text, T_INT, offsetof(CondensedMatrix, data_size), READONLY, PyDoc_STR(bogus_description_text)},
    {NULL}  /* Sentinel */
};

static PyObject* condensedMatrix_get_number_of_rows(CondensedMatrix* self, PyObject *args){
	return Py_BuildValue("i", self->row_length);
}

static PyObject* condensedMatrix_get_data(CondensedMatrix* self, PyObject *args){
	Py_INCREF(self->numpy_array);
	return  PyArray_Return((PyArrayObject*) self->numpy_array);
}

#include "Matrix.Statistics.cpp"

#include "Matrix.Neighbours.cpp"

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
	{"calculateMax", 		(PyCFunction) condensedMatrix_get_max, METH_NOARGS,PyDoc_STR("description")},
	{"calculateMin", 		(PyCFunction) condensedMatrix_get_min, METH_NOARGS,PyDoc_STR("description")},

	// Matrix as graph
	{"get_neighbors_for_node", (PyCFunction)condensedMatrix_get_neighbors_for_node, METH_VARARGS,PyDoc_STR("description")},
	{"choose_node_with_higher_cardinality", (PyCFunction)condensedMatrix_choose_node_with_higher_cardinality, METH_VARARGS,PyDoc_STR("description")},
	{"element_neighbors_within_radius",(PyCFunction)condensedMatrix_get_neighbors_of_node_for_radius, METH_VARARGS,PyDoc_STR("description")},
	//{"calculate_rw_laplacian", (PyCFunction)condensedMatrix_calculate_rw_laplacian, METH_NOARGS,PyDoc_STR("description")},
	//{"calculate_affinity_matrix", (PyCFunction)condensedMatrix_calculate_affinity_matrix, METH_VARARGS,PyDoc_STR("description")},
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
		self->data[pos] = (float) PyFloat_AsDouble(v);
		/////////////////////////////////////////
		// BEWARE!!!!!! SLOW AND REDUNDANT HACK
		/////////////////////////////////////////
		// ALTERNATIVE: in construction use PyObject *PyArray_FROM_OTF(PyObject* obj, int typenum, int requirements) with
		// http://docs.scipy.org/doc/numpy/user/c-info.how-to-extend.html?highlight=pyarray_simplenew#PyArray_SimpleNew
		// NPY_ARRAY_INOUT_ARRAY
		// NPY_ARRAY_FORCECAST
		// NPY_ARRAY_ENSURECOPY
		// NPY_ARRAY_ENSUREARRAY
		if(self->numpy_array != NULL){
			((double*) PyArray_GETPTR1(self->numpy_array,0))[pos] =self->data[pos];
		}
	}
	else{
		if (i!=j){
			pos = calc_vector_pos(j,i,self);
			self->data[pos] = (float) PyFloat_AsDouble(v);
			/////////////////////////////////////////
			// BEWARE!!!!!! SLOW AND REDUNDANT HACK
			/////////////////////////////////////////
			if(self->numpy_array != NULL){
				((double*) PyArray_GETPTR1(self->numpy_array,0))[pos] =self->data[pos];
			}
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
#if PY_VERSION_HEX >= 0x03000000
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                         						/*ob_size*/
#endif
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


#if PY_MAJOR_VERSION >= 3
#define MODULE_INIT_ERROR NULL
PyMODINIT_FUNC PyInit_condensedMatrix(void){
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "condensedMatrix",           /* m_name */
        "Fast Access Condensed Matrix",  /* m_doc */
        -1,                    /* m_size */
        NULL,                  /* m_methods */
        NULL,                  /* m_reload */
        NULL,                  /* m_traverse */
        NULL,                  /* m_clear */
        NULL,                  /* m_free */
    };
#else
#define MODULE_INIT_ERROR
PyMODINIT_FUNC initcondensedMatrix(void){
#endif

    PyObject* module;

    if (PyType_Ready(&CondensedMatrixType) < 0)
        return MODULE_INIT_ERROR;

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&moduledef);
#else
    module = Py_InitModule3("condensedMatrix", NULL,"Fast Access Condensed Matrix");
#endif
    if (module == NULL)
          return MODULE_INIT_ERROR;

    Py_INCREF(&CondensedMatrixType);
    PyModule_AddObject(module, "CondensedMatrix", (PyObject*) &CondensedMatrixType);

    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
