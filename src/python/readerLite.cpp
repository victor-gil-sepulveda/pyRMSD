#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <iostream>
#include <cstring>
#include "../pdbreaderlite/PDBReader.h"
using namespace std;

void parse_params(PyObject *args, char*& path, char*& atom_name_filter){
	int filter_length = 0;
	char* path_c = NULL;
	char* name_filter_c = NULL;
	atom_name_filter = NULL;

	if (!PyArg_ParseTuple(args, "s|s#",&path_c, &name_filter_c, &filter_length)){
    	PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters");
    }

	int path_len = strlen ( path_c );
	path = (char*) malloc(path_len*sizeof(char));
	strcpy ( path, path_c);

	if (filter_length != 0 && filter_length!=4){
		cout<<"Atom name selector had "<<filter_length<<" characters."<<endl;
		PyErr_SetString(PyExc_RuntimeError, "If specified, the atom name must have exactly 4 characters.");
	}
	else{
		if (filter_length != 0){
			atom_name_filter = (char*) malloc(filter_length*sizeof(char));
			strcpy( atom_name_filter, name_filter_c);
		}
	}
}

PyArrayObject* pack_results(vector<double>& coordinates, int number_of_confomations, int number_of_atoms){
	PyArrayObject* coordinates_numpy;
	
	//Store mem for this guys (it will be erased by the python object)
	double * object_coords = new double[coordinates.size()];
	for (unsigned int i =0; i<coordinates.size(); i++){
		object_coords[i] = coordinates[i];
	}

	npy_intp dims[1] = {coordinates.size()};
	coordinates_numpy = (PyArrayObject *) PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, object_coords);

    return coordinates_numpy;
}

static PyObject* readPDB(PyObject *self, PyObject *args){
	char* path = NULL;
	char* atom_name_filter = NULL;

	parse_params(args,path,atom_name_filter);

	PDBReader reader;

	if (!atom_name_filter){
		reader.read(path, "");
	}
	else{
		reader.read(path, atom_name_filter);
	}

	PyArrayObject* coords_array = pack_results(reader.all_coordinates, reader.number_of_models, reader.number_of_atoms);
    return Py_BuildValue("(Oii)", coords_array, reader.number_of_models, reader.number_of_atoms);
}

static PyMethodDef pdbReaderMethods[] = {
    {"readPDB",  readPDB, METH_VARARGS,""},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initpdbReader(void){
    (void) Py_InitModule("pdbReader", pdbReaderMethods);

    import_array();
}
