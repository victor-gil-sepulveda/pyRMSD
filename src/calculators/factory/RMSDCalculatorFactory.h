/*
 * RMSDCalculatorFactory.h
 *
 *  Created on: 06/03/2013
 *      Author: victor
 */

#ifndef RMSDCALCULATORFACTORY_H_
#define RMSDCALCULATORFACTORY_H_
#include "RMSDCalculatorTypes.h"

class RMSDCalculator;

class RMSDCalculatorFactory {
	public:
		RMSDCalculatorFactory();
		virtual ~RMSDCalculatorFactory();

		static RMSDCalculator* createCalculator(
				RMSDCalculatorType type,
				int numberOfConformations,
				int atomsPerConformation,
				double* allCoordinates,
				int atomsPerCalculationConformation = 0,
				double* allCalculationCoordinates = NULL,
				int number_of_threads = 4,
				int threads_per_block = 8,
				int blocks_per_grid = 16,
				bool modify_coordinates = false);
};

#endif /* RMSDCALCULATORFACTORY_H_ */
