/*
 * KABSCHSerialKernel.h
 *
 *  Created on: 14/10/2013
 *      Author: victor
 */

#ifndef NOSUPSERIALKERNEL_H_
#define NOSUPSERIALKERNEL_H_

#include "../KernelFunctions.h"

class NoSuperpositionSerialKernel: public KernelFunctions {
	public:
		NoSuperpositionSerialKernel();
		virtual ~NoSuperpositionSerialKernel();

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
};

#endif /* KABSCHSERIALKERNEL_H_ */
