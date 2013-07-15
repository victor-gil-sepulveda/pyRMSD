/*
 * RMSDCalculationData.h
 *
 *  Created on: Jul 13, 2013
 *      Author: victor
 */

#ifndef RMSDCALCULATIONDATA_H_
#define RMSDCALCULATIONDATA_H_

class RMSDCalculationData {
	public:
		RMSDCalculationData(int numberOfConformations,
							int atomsPerFittingConformation,
							double* fittingCoordinates,
							int atomsPerCalculationConformation,
							double* calculationCoordinates);

		virtual ~RMSDCalculationData();

		bool hasCalculationCoordinatesSet();

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
};
#endif /* RMSDCALCULATIONDATA_H_ */
