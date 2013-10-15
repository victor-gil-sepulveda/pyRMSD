/*
 * KABSCHSerialKernel.h
 *
 *  Created on: 14/10/2013
 *      Author: victor
 */

#ifndef NOSUPOMPKERNEL_H_
#define NOSUPOMPKERNEL_H_

#include "../KernelFunctions.h"

class NoSuperpositionOmpKernel: public KernelFunctions {
	public:
		NoSuperpositionOmpKernel(int number_of_threads);
		virtual ~NoSuperpositionOmpKernel();

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
};

#endif /* KABSCHSERIALKERNEL_H_ */
