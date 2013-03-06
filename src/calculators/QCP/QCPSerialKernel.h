#ifndef _KERNEL_SERIAL_H_
#define _KERNEL_SERIAL_H_

#include "../KernelFunctions.h"
#include <cstdlib>

class ThRMSDSerialKernel: public KernelFunctions{

	public:
		ThRMSDSerialKernel(){}
		virtual ~ThRMSDSerialKernel(){}

		double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords,
				int number_of_atoms, double* rot_matrix = NULL);

		virtual void oneVsFollowingFitEqualCalcWithoutConfRotation(
					double* reference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates);

		virtual void oneVsFollowingFitEqualCalcWithConfRotation(
					double* reference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates);

		virtual void oneVsFollowingFitDiffersCalcWithoutConfRotation(
					double* fitReference,
					double* calcReference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates,
					int coordinatesPerRMSDConformation,
					int atomsPerRMSDConformation,
					double *allRMSDCoordinates);

		virtual void oneVsAllFitDiffersCalcWithConfRotation(
					double* fitReference,
					double* calcReference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates,
					int coordinatesPerRMSDConformation,
					int atomsPerRMSDConformation,
					double *allRMSDCoordinates);

		double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms, double* rot_matrix = NULL);


};
#endif
