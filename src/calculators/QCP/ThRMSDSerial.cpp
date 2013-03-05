#include "ThRMSDSerial.h"
#include "kernel/kernel_functions_serial.h"
#include "../RMSDTools.h"
#include <assert.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
///	Class constructor. Allocates memory and copies coordinates with the required precision type.
///
/// \param 	numberOfConformations [In] Total number of conformations in the coordinates array.
///
/// \param 	atomsPerConformation [In] Number of atoms of every conformation.
///
/// \param 	coords [In] Coordinates array with numberOfConformations*atomsPerConformation*3 elements.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
QCPSerialCalculator::QCPSerialCalculator(int numberOfConformations, int atomsPerConformation, double* coords):
				RMSDCalculator(numberOfConformations, atomsPerConformation, coords){
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class destructor. Frees memory.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
QCPSerialCalculator::~QCPSerialCalculator(){}

KernelFunctions* QCPSerialCalculator::getKernelFunctions(){
	if (this->kernelFunctions == NULL){
		this->kernelFunctions =  new ThRMSDSerialKernel;
	}

	return this->kernelFunctions;
}
