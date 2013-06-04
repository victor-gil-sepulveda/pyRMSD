#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <iostream>
#include "../calculators/factory/RMSDCalculatorFactory.h"
#include "../calculators/RMSDCalculator.h"
using namespace std;


void parse_params_for_one_vs_others(
		PyObject *args,
		RMSDCalculatorType& calculatorType,
		double*& coords_list, int& number_of_atoms,
		double*& calc_coords_list, int& number_of_calc_atoms,
		int& conformation_number, int& number_of_conformations,
    	int& number_of_threads,	int& threads_per_block,	int& blocks_per_grid,
    	int& modify_coordinates){
	
	PyArrayObject* coords_list_obj,
				 * calc_coords_list_obj;

	if (!PyArg_ParseTuple(
			args,
			"iO!iO!iiiiiii",
			&calculatorType,
			&PyArray_Type, &coords_list_obj, &number_of_atoms,
			&PyArray_Type, &calc_coords_list_obj, &number_of_calc_atoms,
			&conformation_number, &number_of_conformations,
			&number_of_threads, &threads_per_block, &blocks_per_grid,
			&modify_coordinates)){

		PyErr_SetString(
				PyExc_RuntimeError,
				"Error parsing parameters"
		);
    }

	if (coords_list_obj->descr->type_num != NPY_DOUBLE || coords_list_obj->nd != 1){
    	PyErr_SetString(PyExc_RuntimeError, "First parameters must be a double array.");
    }

    // Get the data pointer
    coords_list = (double *) (coords_list_obj->data);

    if(number_of_calc_atoms != 0){
    	calc_coords_list = (double *) (calc_coords_list_obj->data);
    }
    else{
    	calc_coords_list = NULL;
    }

    // Length of the vector
    int coords_list_len = coords_list_obj->dimensions[0];

    // Check the length
    if (coords_list_len == 0 || coords_list_len%3 != 0){
		PyErr_SetString(PyExc_RuntimeError, "Invalid size for the coordinates parameter.");
	}
    //TODO check length of calc coordinates
}

void parse_params_for_condensed_matrix(
		PyObject *args,
		RMSDCalculatorType& calculatorType,
		double*& coords_list, int& number_of_atoms,
		double*& calc_coords_list, int& number_of_calc_atoms,
		int& number_of_conformations,
		int& number_of_threads,	int& threads_per_block,	int& blocks_per_grid){

	PyArrayObject* coords_list_obj,
	 	 	 	 * calc_coords_list_obj;

	if (!PyArg_ParseTuple(
			args,
			"iO!iO!iiiii",
			&calculatorType,
			&PyArray_Type, &coords_list_obj, &number_of_atoms,
			&PyArray_Type, &calc_coords_list_obj, &number_of_calc_atoms,
			&number_of_conformations,
			&number_of_threads, &threads_per_block, &blocks_per_grid)){
    	PyErr_SetString(PyExc_RuntimeError, "Error parsing parameters");
    }

    if (coords_list_obj->descr->type_num != NPY_DOUBLE || coords_list_obj->nd != 1){
    	PyErr_SetString(PyExc_RuntimeError, "First parameters must be a double array.");
    }
    // Get the data pointer
    coords_list = (double *) (coords_list_obj->data);
    if(number_of_calc_atoms != 0){
		calc_coords_list = (double *) (calc_coords_list_obj->data);
	}
	else{
		calc_coords_list = NULL;
	}
    coords_list = (double *) (coords_list_obj->data);
    // Length of the vector
    int coords_list_len = coords_list_obj->dimensions[0];

    // Check the length
    if (coords_list_len == 0 || coords_list_len%3 != 0){
		PyErr_SetString(PyExc_RuntimeError, "Invalid size for the coordinates parameter.");
	}
}

PyArrayObject* embed_rmsd_data(vector<double>& rmsd){
	npy_intp dims[1] = {rmsd.size()};
	PyArrayObject* rmsds_list_obj = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

	double* rmsd_data = (double*) (rmsds_list_obj->data);
	for (unsigned int i =0; i<rmsd.size(); i++){
		rmsd_data[i] = rmsd[i];
	}
	
	return rmsds_list_obj;
}

static PyObject* oneVsFollowing(PyObject *self, PyObject *args){

	int conformation_number;
	int atoms_per_conformation;
	int atoms_per_calc_conformation;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;
	int modify_coordinates;

	RMSDCalculatorType calculatorType;
	double* all_coordinates;
	double* all_calc_coordinates;
	vector<double> rmsd;

	parse_params_for_one_vs_others(
			args,
			calculatorType,
			all_coordinates, atoms_per_conformation,
			all_calc_coordinates, atoms_per_calc_conformation,
			conformation_number, number_of_conformations,
			number_of_threads, threads_per_block, blocks_per_grid,
			modify_coordinates);

	rmsd.resize(number_of_conformations-conformation_number-1);

	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
					calculatorType,
					number_of_conformations,
					atoms_per_conformation,
					all_coordinates,
					atoms_per_calc_conformation,
					all_calc_coordinates,
					number_of_threads,
					threads_per_block,
					blocks_per_grid,
					(bool) modify_coordinates);

    rmsdCalculator->oneVsFollowing(conformation_number,&(rmsd[0]));

	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);

	delete rmsdCalculator;

	return PyArray_Return(rmsds_list_obj);
}

static PyObject* iterativeSuperposition(PyObject *self, PyObject *args){
	int atoms_per_conformation;
	int atoms_per_calc_conformation;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;
	double* all_coordinates;
	double* all_calc_coordinates;
	vector<double> rmsd;
	RMSDCalculatorType calculatorType;


	parse_params_for_condensed_matrix(
			args,
			calculatorType,
			all_coordinates, atoms_per_conformation,
			all_calc_coordinates, atoms_per_calc_conformation,
			number_of_conformations,
			number_of_threads, threads_per_block,blocks_per_grid
	);

	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
				calculatorType,
				number_of_conformations,
				atoms_per_conformation,
				all_coordinates,
				atoms_per_calc_conformation,
				all_calc_coordinates,
				number_of_threads,
				threads_per_block,
				blocks_per_grid);

	rmsdCalculator->iterativeSuperposition(1e-4);

	delete rmsdCalculator;

	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* calculateRMSDCondensedMatrix(PyObject *self, PyObject *args){
	int atoms_per_conformation;
	int atoms_per_calc_conformation;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;
	double* all_coordinates;
	double* all_calc_coordinates;
	vector<double> rmsd;
	RMSDCalculatorType calculatorType;


	parse_params_for_condensed_matrix(
			args,
			calculatorType,
			all_coordinates, atoms_per_conformation,
			all_calc_coordinates, atoms_per_calc_conformation,
			number_of_conformations,
			number_of_threads, threads_per_block,blocks_per_grid
	);

	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
			calculatorType,
			number_of_conformations,
			atoms_per_conformation,
			all_coordinates,
			atoms_per_calc_conformation,
			all_calc_coordinates,
			number_of_threads,
			threads_per_block,
			blocks_per_grid);

	rmsdCalculator->calculateRMSDCondensedMatrix(rmsd);

	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);
	delete rmsdCalculator;

	return PyArray_Return(rmsds_list_obj);
}

static PyMethodDef pyRMSDMethods[] = {
    {"oneVsFollowing",  oneVsFollowing, METH_VARARGS,""},
    {"calculateRMSDCondensedMatrix",  calculateRMSDCondensedMatrix, METH_VARARGS,""},
    {"iterativeSuperposition",  iterativeSuperposition, METH_VARARGS,""},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initcalculators(void){
    (void) Py_InitModule("calculators", pyRMSDMethods);

    import_array();
}
