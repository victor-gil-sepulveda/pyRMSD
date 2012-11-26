#include <python2.7/Python.h>
#include <numpy/arrayobject.h>
#include "NumpyHelperFuncs.h"
#include "../serial/RMSDSerial.h"
#include "../serial/RMSDomp.h"
#include "../theobald/ThRMSDCuda.cuh"
#include "../theobald/ThRMSDSerial.h"
#include "../theobald/ThRMSDSerialOmp.h"
#include <omp.h>
#include <vector>
#include <iostream>
using namespace std;

enum calcType{
	SERIAL_CALCULATOR = 0,
	OPENMP_CALCULATOR,
	THEOBALD_CUDA_CALCULATOR,
	THEOBALD_SERIAL_CALCULATOR,
	THEOBALD_SERIAL_OMP_CALCULATOR
};

void parse_params_for_one_vs_others(PyObject *args, calcType& cType,
    		int& number_of_atoms, int& conformation_number, int& number_of_conformations,
    		int& number_of_threads, int& threads_per_block,int& blocks_per_grid,
    		double*& coords_list){
	
	PyArrayObject* coords_list_obj;
	
	if (!PyArg_ParseTuple(args, "iO!iiiiii",&cType,&PyArray_Type, &coords_list_obj,
	&number_of_atoms, &conformation_number, &number_of_conformations,
	&number_of_threads, &threads_per_block, &blocks_per_grid)){
    	PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters");
    }

    if (not_doublevector(coords_list_obj)){
    	PyErr_SetString(PyExc_RuntimeError, "First parameters must be a double array.");
    }

    // Get the data pointer
    coords_list = pyvector_to_Carrayptrs(coords_list_obj);

    // Length of the vector
    int coords_list_len = coords_list_obj->dimensions[0];

    // Check the length
    if (coords_list_len == 0 || coords_list_len%3 != 0){
		PyErr_SetString(PyExc_RuntimeError, "Invalid size for the coordinates parameter.");
	}
}

void parse_params_for_condensed_matrix(PyObject *args,calcType& cType,
    		int& number_of_atoms, int& number_of_conformations,
    		int& number_of_threads, int& threads_per_block,int& blocks_per_grid,
    		double*& coords_list){

	PyArrayObject* coords_list_obj;

	if (!PyArg_ParseTuple(args, "iO!iiiii",&cType,&PyArray_Type, &coords_list_obj,
	&number_of_atoms, &number_of_conformations,
	&number_of_threads, &threads_per_block, &blocks_per_grid)){
    	PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters");
    }

    if (not_doublevector(coords_list_obj)){
    	PyErr_SetString(PyExc_RuntimeError, "First parameters must be a double array.");
    }

    // Get the data pointer
    coords_list = pyvector_to_Carrayptrs(coords_list_obj);

    // Length of the vector
    int coords_list_len = coords_list_obj->dimensions[0];

    // Check the length
    if (coords_list_len == 0 || coords_list_len%3 != 0){
		PyErr_SetString(PyExc_RuntimeError, "Invalid size for the coordinates parameter.");
	}
}

PyArrayObject* embed_rmsd_data(vector<double>& rmsd){
	PyArrayObject* rmsds_list_obj;
	// Possible leak?
	double* rmsd_data = new double[rmsd.size()];
	for (unsigned int i =0; i<rmsd.size(); i++){
		rmsd_data[i] = rmsd[i];
	}
	
	npy_intp dims[1] = {rmsd.size()};
    rmsds_list_obj = (PyArrayObject *) PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,rmsd_data);
	return rmsds_list_obj;
}

RMSD* getCalculator(calcType cType, int numberOfConformations, int atomsPerConformation,
		int number_of_threads, int threads_per_block, int blocks_per_grid,
		double* Coordinates){

	switch(cType){

		case SERIAL_CALCULATOR:
			return new RMSDSerial(numberOfConformations,atomsPerConformation,Coordinates);
			break;

		case OPENMP_CALCULATOR:
			return new RMSDomp(numberOfConformations,atomsPerConformation,Coordinates);
			break;

#ifdef USE_CUDA
		case THEOBALD_CUDA_CALCULATOR:
			return new ThRMSDCuda(numberOfConformations, atomsPerConformation, Coordinates, threads_per_block, blocks_per_grid);
			break;
#endif
		case THEOBALD_SERIAL_CALCULATOR:
			return new ThRMSDSerial(numberOfConformations,atomsPerConformation,Coordinates);
			break;

		case THEOBALD_SERIAL_OMP_CALCULATOR:
			return new ThRMSDSerialOmp(numberOfConformations,atomsPerConformation,Coordinates);
			break;

		default:
			return NULL;
	}
}

static PyObject* oneVsFollowing(PyObject *self, PyObject *args){

	int conformation_number;
	int atoms_per_conformation;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;

	calcType cType;
	double* all_coordinates;
	vector<double> rmsd;

	parse_params_for_one_vs_others(args, cType, atoms_per_conformation, conformation_number, number_of_conformations,
			number_of_threads, threads_per_block, blocks_per_grid, all_coordinates);

	rmsd.resize(number_of_conformations-conformation_number-1);

    RMSD* rmsdCalculator = getCalculator(cType, number_of_conformations, atoms_per_conformation,
    		number_of_threads, threads_per_block, blocks_per_grid, all_coordinates);

    rmsdCalculator->oneVsFollowing(conformation_number,&(rmsd[0]));

	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);
	delete rmsdCalculator;

    return PyArray_Return(rmsds_list_obj);
}

static PyObject* calculateRMSDCondensedMatrix(PyObject *self, PyObject *args){
	int atoms_per_conformation;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;
	double* all_coordinates;
	calcType cType;
	vector<double> rmsd;

	parse_params_for_condensed_matrix(args, cType, atoms_per_conformation, number_of_conformations,
			number_of_threads, threads_per_block, blocks_per_grid, all_coordinates);

	RMSD* rmsdCalculator = getCalculator(cType, number_of_conformations, atoms_per_conformation,
			number_of_threads, threads_per_block, blocks_per_grid, all_coordinates);

	rmsdCalculator->calculateRMSDCondensedMatrix(rmsd);

	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);
	delete rmsdCalculator;

	return PyArray_Return(rmsds_list_obj);
}

static PyMethodDef pyRMSDMethods[] = {
    {"oneVsFollowing",  oneVsFollowing, METH_VARARGS,""},
    {"calculateRMSDCondensedMatrix",  calculateRMSDCondensedMatrix, METH_VARARGS,""},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcalculators(void){
    (void) Py_InitModule("calculators", pyRMSDMethods);

    import_array();
}
