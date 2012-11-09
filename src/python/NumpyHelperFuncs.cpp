// FROM http://www.scipy.org/Cookbook/C_Extensions/NumPy_arrays
#include "NumpyHelperFuncs.h"
#include <python2.7/Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
return 1 if an error and raise exception */
int  not_doublevector(PyArrayObject *vec){
	if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {
		PyErr_SetString(PyExc_ValueError,
				"In not_doublevector: array must be of type Double and 1 dimensional (n).");
		return 1;
	}
	return 0;
}


/* ==== Create 1D Carray from PyArray ======================
     Assumes PyArray is contiguous in memory.             */
double *pyvector_to_Carrayptrs(PyArrayObject *arrayin){
	return (double *) arrayin->data;  /* pointer to arrayin data as double */
}

