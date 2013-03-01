#include "ThRMSDSerialOmp.h"
#include "ThRMSDSerial.h"
#include "kernel_functions_omp.h"
#include <assert.h>
#include <iostream>
#include <omp.h>

using namespace std;

ThRMSDSerialOmp::ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation,
		double* coords, int omp_threads):
				ThRMSDSerial(numberOfConformations, atomsPerConformation, coords){

	this->omp_threads = omp_threads;
	omp_set_num_threads(this->omp_threads);

	this->kernelFunctions = new ThRMSDSerialOmpKernel;
}

ThRMSDSerialOmp::~ThRMSDSerialOmp(){}
