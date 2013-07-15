/*
 * RMSDCalculationData.cpp
 *
 *  Created on: Jul 13, 2013
 *      Author: victor
 */

#include "RMSDCalculationData.h"
#include <stdlib.h>

RMSDCalculationData::RMSDCalculationData(int numberOfConformations,
											int atomsPerFittingConformation,
											double* fittingCoordinates,
											int atomsPerCalculationConformation,
											double* calculationCoordinates){

	this->numberOfConformations = numberOfConformations;

	this->atomsPerFittingConformation = atomsPerFittingConformation;
	this->fittingConformationLength = atomsPerFittingConformation*3;
	this->fittingCoordinatesLength = atomsPerFittingConformation*numberOfConformations*3;
	this->fittingCoordinates = fittingCoordinates;

	this->atomsPerCalculationConformation = atomsPerCalculationConformation;
	this->calculationConformationLength = atomsPerCalculationConformation*3;
	this->calculationCoordinatesLength = atomsPerCalculationConformation*numberOfConformations*3;
	this->calculationCoordinates = calculationCoordinates;
}

RMSDCalculationData::~RMSDCalculationData(){}

bool RMSDCalculationData::hasCalculationCoordinatesSet(){
	return this->calculationCoordinates != NULL;
}
