#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include "PDBReader.h"
using namespace std;

#define calc_vector_pos(i,j,matrix) (i*(matrix->row_length_minus_one) - i - ((( i-1)*i) >> 1 ) + j - 1)

typedef struct {
    PyObject_HEAD
    PDBReader* reader;
//    double* coordinates;
} pdbreader;


static PyMemberDef pdbreader_members[] = {
    {NULL}  /* Sentinel */
};

static PyObject* pdbreader_new(PyTypeObject *type, PyObject *args, PyObject *kwds){
	pdbreader* self = (pdbreader*) type->tp_alloc(type, 0);
	if(self != NULL){
		self->reader = new PDBReader;
//		self->coordinates = NULL;
	}

    return (PyObject *) self;
}

static void pdbreader_dealloc(pdbreader* self){
	if(self->reader!=NULL){
		delete self->reader;
	}

//	if(self->coordinates != NULL){
//		delete [] self->coordinates;
//	}
}

static int pdbreader_init(pdbreader *self, PyObject *args, PyObject *kwds){
	return 0;
}

void parse_params(PyObject *args, char*& path, char*& atom_name_filter){
	int filter_length = 0;
	char* name_filter_c = NULL;
	atom_name_filter = NULL;

	if (!PyArg_ParseTuple(args, "s|s#",&path, &name_filter_c, &filter_length)){
    	PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters");
    }

	if (filter_length != 0 && filter_length!=4){
		cout<<"Atom name selector had "<<filter_length<<" characters."<<endl;
		PyErr_SetString(PyExc_RuntimeError, "If specified, the atom name must have exactly 4 characters.");
	}
	else{
		if (filter_length != 0){
			atom_name_filter = new char[filter_length+1];
			strcpy( atom_name_filter, name_filter_c);
		}
	}
}

static PyObject* pdbreader_read(pdbreader* self, PyObject *args){
	char* path = NULL;
	char* atom_name_filter = NULL;

	parse_params(args, path, atom_name_filter);

	if (!atom_name_filter){
		self->reader->read(path, "");
	}
	else{
		self->reader->read(path, atom_name_filter);
		delete [] atom_name_filter;
	}

//	if(self->coordinates != NULL){
//		delete [] self->coordinates;
//	}

	npy_intp dims[] = {self->reader->all_coordinates.size()};
//	self->coordinates = new double[dims[0]];
//	double* remote_coords = &(self->reader->all_coordinates[0]);
//	copy(remote_coords,remote_coords+dims[0],self->coordinates);
	PyObject* coords_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, &(self->reader->all_coordinates[0]));
	return Py_BuildValue("(Oii)", coords_array, self->reader->number_of_models, self->reader->number_of_atoms);

}


static PyMethodDef pdbreader_methods[] = {
	{"read", (PyCFunction) pdbreader_read, METH_VARARGS, PyDoc_STR("description")},
	{NULL}  /* Sentinel */
};

static PyTypeObject pdbreaderType = {
#if PY_VERSION_HEX >= 0x03000000
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                         						/*ob_size*/
#endif
    "PDBReader.Reader",      	/*tp_name*/
    sizeof(pdbreader), 	/*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)pdbreader_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0, 				           /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT| Py_TPFLAGS_BASETYPE , /*tp_flags*/
    "Simple PDB reader",           /* tp_doc */
	0,		               					  /* tp_traverse */
	0,		               					  /* tp_clear */
	0,		               					  /* tp_richcompare */
	0,		              					  /* tp_weaklistoffset */
	0,		               	   /* tp_iter */
	0,		               	   /* tp_iternext */
	pdbreader_methods,   /* tp_methods */
	pdbreader_members,   /* tp_members */
	0,                         /* tp_getset */
	0,                         /* tp_base */
	0,                         /* tp_dict */
	0,                         /* tp_descr_get */
	0,                         /* tp_descr_set */
	0,                         /* tp_dictoffset */
	(initproc)pdbreader_init, /* tp_init */
	0,                         		/* tp_alloc */
	pdbreader_new,        		/* tp_new */
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif


#if PY_MAJOR_VERSION >= 3
#define MODULE_INIT_ERROR NULL
PyMODINIT_FUNC PyInit_pdbReader(void){
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pdbReader",           /* m_name */
        "Simple pdb reading",  /* m_doc */
        -1,                    /* m_size */
        NULL,                  /* m_methods */
        NULL,                  /* m_reload */
        NULL,                  /* m_traverse */
        NULL,                  /* m_clear */
        NULL,                  /* m_free */
    };
#else
#define MODULE_INIT_ERROR
PyMODINIT_FUNC initpdbReader(void){
#endif

    PyObject* module;

    if (PyType_Ready(&pdbreaderType) < 0)
        return MODULE_INIT_ERROR;

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&moduledef);
#else
    module = Py_InitModule3("pdbReader", NULL,"Simple pdb reading");
#endif
    if (module == NULL)
          return MODULE_INIT_ERROR;

    Py_INCREF(&pdbreaderType);
    PyModule_AddObject(module, "PDBReader", (PyObject*) &pdbreaderType);

    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
