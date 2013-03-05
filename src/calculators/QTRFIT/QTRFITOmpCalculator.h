/*
 * RMSDomp.h
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#ifndef RMSDOMP_H_
#define RMSDOMP_H_
#include "../RMSDCalculator.h"


class QTRFITOmpCalculator: public RMSDCalculator{

	public:
		QTRFITOmpCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		virtual ~QTRFITOmpCalculator();

	private:
		KernelFunctions* getKernelFunctions();
};

#endif /* RMSDOMP_H_ */
