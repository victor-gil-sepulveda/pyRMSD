/*
 * QTRFITOmpKernel.h
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#ifndef QTRFITOMPKERNEL_H_
#define QTRFITOMPKERNEL_H_

#include "QTRFITSerialKernel.h"

class QTRFITOmpKernel: public QTRFITSerialKernel {
	public:
		QTRFITOmpKernel(int number_of_threads);
		virtual ~QTRFITOmpKernel();

		virtual void oneVsFollowingFitEqualCalcCoords(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				RMSDCalculationData* data);

		virtual void oneVsFollowingFitDiffersCalcCoords(
				double* fitReference,
				double* calcReference,
				int reference_conformation_number,
				double* rmsd,
				RMSDCalculationData* data);

		int number_of_threads;

	private:
		QTRFITOmpKernel(){number_of_threads = 0;};
};

#endif /* QTRFITOMPKERNEL_H_ */
