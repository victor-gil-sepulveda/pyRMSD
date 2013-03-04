#include "ThRMSDSerialOmp.h"
#include "ThRMSDSerial.h"
#include "kernel/kernel_functions_omp.h"
#include <assert.h>
#include <iostream>
#include <omp.h>

using namespace std;

QCPOmpCalculator::QCPOmpCalculator(int numberOfConformations, int atomsPerConformation,
		double* coords, int omp_threads):
				QCPSerialCalculator(numberOfConformations, atomsPerConformation, coords){

	this->omp_threads = omp_threads;
	omp_set_num_threads(this->omp_threads);
}

QCPOmpCalculator::~QCPOmpCalculator(){}

KernelFunctions* QCPOmpCalculator::getKernelFunctions(){
	if (this->kernelFunctions == NULL){
		this->kernelFunctions =  new ThRMSDSerialKernel;
	}

	return this->kernelFunctions;
}
