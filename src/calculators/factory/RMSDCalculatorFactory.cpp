/*
 * RMSDCalculatorFactory.cpp
 *
 *  Created on: 06/03/2013
 *      Author: victor
 */

#include "RMSDCalculatorFactory.h"
#include "../KernelFunctions.h"
#include "../QTRFIT/QTRFITOmpKernel.h"
#include "../QCP/QCPSerialKernel.h"
#include "../QCP/QCPOmpKernel.h"
#include "../RMSDCalculator.h"

#include <iostream>
#include <cstdlib>
#include "../QTRFIT/QTRFITSerialKernel.h"
using namespace std;

RMSDCalculatorFactory::RMSDCalculatorFactory() {}

RMSDCalculatorFactory::~RMSDCalculatorFactory() {}

RMSDCalculator* RMSDCalculatorFactory::createCalculator(
		RMSDCalculatorType type,
		int numberOfConformations,
		int atomsPerConformation,
		double* allCoordinates) {

	KernelFunctions* kernelFunctions;

	switch (type) {
		case KABSCH_SERIAL_CALCULATOR:
					kernelFunctions = NULL;
					break;

		case KABSCH_OMP_CALCULATOR:
					kernelFunctions = NULL;
					break;

		case KABSCH_CUDA_CALCULATOR:
					kernelFunctions = NULL;
					break;

		case QTRFIT_SERIAL_CALCULATOR:
					kernelFunctions = new QTRFITSerialKernel;
					break;

		case QTRFIT_OMP_CALCULATOR:
					kernelFunctions = new QTRFITOmpKernel;
					break;

		case QTRFIT_CUDA_CALCULATOR:
					kernelFunctions = NULL;
					break;

		case QCP_SERIAL_CALCULATOR:
					kernelFunctions = new ThRMSDSerialKernel;
					break;

		case QCP_OMP_CALCULATOR:
					kernelFunctions = new ThRMSDSerialOmpKernel;
					break;

		case QCP_CUDA_CALCULATOR:
					kernelFunctions = NULL;
					break;

		default:
			cout<<"[ERROR] Not kernel type implementation for type: "<<calculatorTypeToString(type)<<endl;
			exit(-1);
	}

	return new RMSDCalculator(numberOfConformations,
			atomsPerConformation, allCoordinates, kernelFunctions);
}
