/*
 * KABSCHSerialKernel.h
 *
 *  Created on: 07/03/2013
 *      Author: victor
 */

#ifndef KABSCHSERIALKERNEL_H_
#define KABSCHSERIALKERNEL_H_

#include "../KernelFunctions.h"

class KABSCHSerialKernel: public KernelFunctions {
	public:
		KABSCHSerialKernel();
		virtual ~KABSCHSerialKernel();

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

		virtual double calculate_rotation_rmsd(
				const double* const referenceCoords,
				const double* const fitCoords,
				int numberOfAtoms,
				double (*const U)[3]);

	protected:
		virtual void calc_correlation_matrix_and_E0(
				double (*const R)[3],
				double* const _E0,
				const double* const referenceCoords,
				const double* const fitCoords,
				int numberOfAtoms);

		virtual bool calculate_rotation_matrix(
				const double (*const R)[3],
				double (*const U)[3],
				double E0,
				double* residual);
};

#endif /* KABSCHSERIALKERNEL_H_ */
