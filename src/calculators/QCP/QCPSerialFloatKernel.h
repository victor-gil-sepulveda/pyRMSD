/*
 * QCPSerialFloatKernel.h
 *
 *  Created on: 18/07/2013
 *      Author: victor
 */

#ifndef QCPSERIALFLOATKERNEL_H_
#define QCPSERIALFLOATKERNEL_H_
#include <stdlib.h>
#include "../KernelFunctions.h"

class QCPSerialFloatKernel: public KernelFunctions {
	public:
		QCPSerialFloatKernel(){};
		virtual ~QCPSerialFloatKernel(){};

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

		float calcRMSDOfTwoConformations( float* first_conformation_coords,
												float* second_conformation_coords,
												int number_of_atoms,
												float* rot_matrix = NULL);

		float innerProduct(float* A,
								float* first_conformation_coords,
								float* second_conformation_coords,
								int number_of_atoms);

		float calcRMSDForTwoConformationsWithTheobaldMethod(float *A,
																	float E0,
																	int number_of_atoms,
																	float* rot_matrix = NULL);
};

#endif /* QCPSERIALFLOATKERNEL_H_ */
