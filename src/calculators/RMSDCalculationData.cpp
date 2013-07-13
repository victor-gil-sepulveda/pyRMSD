/*
 * RMSDCalculationData.cpp
 *
 *  Created on: Jul 13, 2013
 *      Author: victor
 */

#include "RMSDCalculationData.h"

RMSDCalculationData::RMSDCalculationData(int numberOfConformations,
											int atomsPerFittingConformation,
											double* const fittingCoordinates,
											int atomsPerCalculationConformation,
											double* const calculationCoordinates){

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

RMSDCalculationData::~RMSDCalculationData() {}

