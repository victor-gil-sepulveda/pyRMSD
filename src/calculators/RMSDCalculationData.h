/*
 * RMSDCalculationData.h
 *
 *  Created on: Jul 13, 2013
 *      Author: victor
 */

#ifndef RMSDCALCULATIONDATA_H_
#define RMSDCALCULATIONDATA_H_

#include "symmGroups.h"
#include <cstddef>

class RMSDCalculationData {
	public:
		RMSDCalculationData(int numberOfConformations,
							int atomsPerFittingConformation,
							double* fittingCoordinates,
							int atomsPerCalculationConformation,
							double* calculationCoordinates,
							symmGroups* symmetryGroups);

		virtual ~RMSDCalculationData();

		inline bool hasCalculationCoordinatesSet(){
			return this->calculationCoordinates != NULL;
		}

		inline bool hasSymmetryGroups(){
			return this->symmetryGroups != NULL && this->symmetryGroups->size() > 0;
		}

		inline double* getFittingConformationAt(int index){
			return &(this->fittingCoordinates[index*this->fittingConformationLength]);
		}

		inline double* getCalculationConformationAt(int index){
			return &(this->calculationCoordinates[index*this->calculationConformationLength]);
		}

		int numberOfConformations;

		int atomsPerFittingConformation;
		int fittingConformationLength;
		int fittingCoordinatesLength;
		double* fittingCoordinates;

		int atomsPerCalculationConformation;
		int calculationConformationLength;
		int calculationCoordinatesLength;
		double* calculationCoordinates;

		symmGroups* symmetryGroups;
};
#endif /* RMSDCALCULATIONDATA_H_ */
