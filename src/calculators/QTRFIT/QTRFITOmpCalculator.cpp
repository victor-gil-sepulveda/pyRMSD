/*
 * RMSDomp.cpp
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#include "QTRFITOmpCalculator.h"
#include "kernel/QTRFITSerialKernel.h"


QTRFITOmpCalculator::QTRFITOmpCalculator(int numberOfConformations, int atomsPerConformation, double *allCoordinates):
		RMSDCalculator(numberOfConformations, atomsPerConformation,allCoordinates){
}

QTRFITOmpCalculator::~QTRFITOmpCalculator(){}

KernelFunctions* QTRFITOmpCalculator::getKernelFunctions(){
	if(this->kernelFunctions == NULL){
			this->kernelFunctions =  new QTRFITSerialKernel;
	}
	return this->kernelFunctions;
}

