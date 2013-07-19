#ifndef _KERNEL_SERIAL_H_
#define _KERNEL_SERIAL_H_

#include "../KernelFunctions.h"
#include <cstdlib>

class QCPSerialKernel: public KernelFunctions{

	public:
		QCPSerialKernel(){}
		virtual ~QCPSerialKernel(){}


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

		double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords,
				int number_of_atoms, double* rot_matrix = NULL);

		double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms, double* rot_matrix = NULL);


};
#endif
