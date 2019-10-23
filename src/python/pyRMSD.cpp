#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <iostream>
#include "../calculators/factory/RMSDCalculatorFactory.h"
#include "../calculators/RMSDCalculator.h"
#include "../calculators/symmGroups.h"

#if PY_VERSION_HEX >= 0x03000000
#define PyInt_AsLong(x) PyLong_AsLong(x)
#endif

using namespace std;

struct Params{
	RMSDCalculatorType calculator_type;
	double* fitting_coords_list;
	int number_of_fitting_atoms;
	double* calc_coords_list;
	int number_of_calc_atoms;
	int conformation_number;
	int number_of_conformations;
	int number_of_threads;
	int threads_per_block;
	int blocks_per_grid;
	symmGroups* symmetry_groups;
};


PyObject* get_item(PyObject* list_or_tuple, int item){
	
	if(PyTuple_Check(list_or_tuple)){
		// Then is a tuple (and item must be 0 or 1
		return PyTuple_GetItem((PyObject*) list_or_tuple, item);
	}
	else{
		// It is a list
		return PyList_GetItem((PyObject*) list_or_tuple, item);
	}
}
/**
 * Parses a 'symmetry groups' structure to recreate its homologous C structure using vectors.
 *
 * \param list_obj Is the python structure (python array of pairs (also tuples) of n-ary tuples.
 *
 * \param symmetry_groups C conversion of this structure using STL objects.
 */
void parse_symmetry_groups(PyListObject* list_obj, symmGroups& symmetry_groups){
	// Precondition: symmetry_groups is an empty vector
	// and symmetry groups have the correct structure and size
	
	int number_of_symmetry_groups = PyList_Size((PyObject*) list_obj);
	//cout<<"DBG: Num. symmgroups "<<number_of_symmetry_groups<<endl;
	
	for (int i = 0; i < number_of_symmetry_groups; ++i){
		PyObject* symm_group =  get_item((PyObject*)list_obj, i);
		vector<pair<int, int> > symm_group_vector;
		int symmetry_elements = PyList_Size( symm_group);
		for (int j=0; j < symmetry_elements; ++j){
			PyObject* symm_pair = get_item(symm_group, j);
			symm_group_vector.push_back(pair<int, int> (
														(int) PyInt_AsLong(get_item(symm_pair, 0)),
														(int) PyInt_AsLong(get_item(symm_pair, 1))
														)
														);
		}
	
		symmetry_groups.push_back(symm_group_vector);
	}
}

/**
 * Parses the parameters given to any of the Python functions and fills a C parameters structure.
 *
 * \params args Python tuple containing function arguments.
 *
 * \params params C structure which will store the parameters.
 *
 */
void parse_params(PyObject *args, Params* params){

	PyArrayObject *fit_coords_list_obj,
				  *calc_coords_list_obj;

	PyListObject* symmetry_groups_list_obj;

	// We can discriminate the type of function by the number of parameters
	int num_params = PyTuple_Size(args);
	bool parsing_ok = false;

	switch (num_params){
		case 11: // one vs following
			parsing_ok = (bool) PyArg_ParseTuple(		args,
														"iO!iO!iiiO!iii",
														&(params->calculator_type),
														&PyArray_Type,
														&fit_coords_list_obj,
														&(params->number_of_fitting_atoms),
														&PyArray_Type,
														&calc_coords_list_obj,
														&(params->number_of_calc_atoms),
														&(params->conformation_number),
														&(params->number_of_conformations),
														&PyList_Type,
														&symmetry_groups_list_obj,
														&(params->number_of_threads),
														&(params->threads_per_block),
														&(params->blocks_per_grid));
			break;
		case 10: // matrix and iterative superposition
			parsing_ok = (bool) PyArg_ParseTuple(		args,
														"iO!iO!iiO!iii",
														&(params->calculator_type),
														&PyArray_Type,
														&fit_coords_list_obj,
														&(params->number_of_fitting_atoms),
														&PyArray_Type,
														&calc_coords_list_obj,
														&(params->number_of_calc_atoms),
														&(params->number_of_conformations),
														&PyList_Type,
														&symmetry_groups_list_obj,
														&(params->number_of_threads),
														&(params->threads_per_block),
														&(params->blocks_per_grid));
			break;
		default:
			PyErr_SetString(PyExc_RuntimeError,
							"Unexpected number of parameters.");
	}

	if (!parsing_ok){
		PyErr_SetString(PyExc_RuntimeError,
						"Error parsing parameters.");
	}

	if (fit_coords_list_obj->descr->type_num != NPY_DOUBLE || fit_coords_list_obj->nd != 1){
		PyErr_SetString(PyExc_RuntimeError,
						"First parameters must be a double array.");
	}

	// Get the data pointer
	params->fitting_coords_list = (double *) (fit_coords_list_obj->data);

	if(params->number_of_calc_atoms != 0){
		params->calc_coords_list = (double *) (calc_coords_list_obj->data);
	}
	else{
		params->calc_coords_list = NULL;
	}

	// Parsing symmetry groups
	params->symmetry_groups = new symmGroups;
	parse_symmetry_groups(symmetry_groups_list_obj, *(params->symmetry_groups));
}

PyArrayObject* embed_rmsd_data(vector<double>& rmsd){
	npy_intp dims[1] = {rmsd.size()};
	PyArrayObject* rmsds_list_obj = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

	double* rmsd_data = (double*) (rmsds_list_obj->data);
	for (unsigned int i = 0; i < rmsd.size(); ++i){
		rmsd_data[i] = rmsd[i];
	}
	
	return rmsds_list_obj;
}

static PyObject* oneVsFollowing(PyObject *self, PyObject *args){

	Params params;
	vector<double> rmsd;

	parse_params(args, &params);

	rmsd.resize(params.number_of_conformations-params.conformation_number-1);
	
	//cout<<"DBG: Parsed_size "<<params.symmetry_groups->size()<<endl;
	//for (unsigned i=0;i<params.symmetry_groups->size();++i){
	//	vector<pair<int,int> > group = params.symmetry_groups->at(i);
	//	for (unsigned j=0;j<group.size();++j){
	//		pair<int,int> pair = group[j];
	//		cout<<" ("<<pair.first<<","<<pair.second<<") ";
	//	}
	//}
	//cout<<endl;
	
	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
									params.calculator_type,
									params.number_of_conformations,
									params.number_of_fitting_atoms,
									params.fitting_coords_list,
									params.number_of_calc_atoms,
									params.calc_coords_list,
									params.symmetry_groups,
									params.number_of_threads,
									params.threads_per_block,
									params.blocks_per_grid);

	rmsdCalculator->oneVsFollowing(params.conformation_number,&(rmsd[0]));

	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);

	delete rmsdCalculator;
	delete params.symmetry_groups;

	return PyArray_Return(rmsds_list_obj);
}

static PyObject* iterativeSuperposition(PyObject *self, PyObject *args){
	Params params;
	vector<double> rmsd;

	parse_params(args, &params);

	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
									params.calculator_type,
									params.number_of_conformations,
									params.number_of_fitting_atoms,
									params.fitting_coords_list,
									params.number_of_calc_atoms,
									params.calc_coords_list,
									params.symmetry_groups,
									params.number_of_threads,
									params.threads_per_block,
									params.blocks_per_grid);

	rmsdCalculator->iterativeSuperposition(1e-4);

	delete rmsdCalculator;
	delete params.symmetry_groups;

	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* calculateRMSDCondensedMatrix(PyObject *self, PyObject *args){
	Params params;
	vector<double> rmsd;

	parse_params(args, &params);

//  cout<<"DBG: Parsed_size "<<params.symmetry_groups->size()<<endl;
//	for (unsigned i=0;i<params.symmetry_groups->size();++i){
//		vector<pair<int,int> > group = params.symmetry_groups->at(i);
//		for (unsigned j=0;j<group.size();++j){
//			pair<int,int> pair = group[j];
//			cout<<" ("<<pair.first<<","<<pair.second<<") ";
//		}
//	}
//	cout<<endl;
//	
//	cout<<"DBG: Inner "<<params.number_of_threads<< " "<<params.threads_per_block<<" "<<params.blocks_per_grid<<endl;

	RMSDCalculator* rmsdCalculator = RMSDCalculatorFactory::createCalculator(
												params.calculator_type,
												params.number_of_conformations,
												params.number_of_fitting_atoms,
												params.fitting_coords_list,
												params.number_of_calc_atoms,
												params.calc_coords_list,
												params.symmetry_groups,
												params.number_of_threads,
												params.threads_per_block,
												params.blocks_per_grid);
    
	rmsdCalculator->calculateRMSDCondensedMatrix(rmsd);
    
	PyArrayObject* rmsds_list_obj = embed_rmsd_data(rmsd);

    
	delete rmsdCalculator;
	delete params.symmetry_groups;
	
    //return PyArray_Return(rmsds_list_obj);
    //Py_INCREF(rmsds_list_obj);
	return (PyObject*) rmsds_list_obj;
}

static PyMethodDef pyRMSDMethods[] = {
    {"oneVsFollowing",  oneVsFollowing, METH_VARARGS,""},
    {"calculateRMSDCondensedMatrix",  calculateRMSDCondensedMatrix, METH_VARARGS,""},
    {"iterativeSuperposition",  iterativeSuperposition, METH_VARARGS,""},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_calculators(void){
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "calculators",         /* m_name */
        "calculators",         /* m_doc */
        -1,                    /* m_size */
        pyRMSDMethods,         /* m_methods */
        NULL,                  /* m_reload */
        NULL,                  /* m_traverse */
        NULL,                  /* m_clear */
        NULL,                  /* m_free */
    };
    PyObject *module;
#else
PyMODINIT_FUNC initcalculators(void){
#endif

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&moduledef);
#else
        (void) Py_InitModule("calculators", pyRMSDMethods);
#endif

    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
