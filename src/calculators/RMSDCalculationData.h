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
							double* const fittingCoordinates,
							int atomsPerCalculationConformation,
							double* const calculationCoordinates);

		virtual ~RMSDCalculationData();

		int numberOfConformations;

		int atomsPerFittingConformation;
		int fittingConformationLength;
		int fittingCoordinatesLength;
		double* const fittingCoordinates;

		int atomsPerCalculationConformation;
		int calculationConformationLength;
		int calculationCoordinatesLength;
		double* const calculationCoordinates;
};
#endif /* RMSDCALCULATIONDATA_H_ */
