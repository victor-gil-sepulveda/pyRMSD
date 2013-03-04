#ifndef _KERNEL_SERIAL_H_
#define _KERNEL_SERIAL_H_

#include "KernelFunctions.h"

class ThRMSDSerialKernel: public KernelFunctions{
	public:
		ThRMSDSerialKernel(){}
		virtual ~ThRMSDSerialKernel(){}

		double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms, double* rot_matrix = NULL);

		double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords,
				int number_of_atoms, double* rot_matrix = NULL);

		virtual void calcRMSDOfOneVsFollowing(double* all_coordinates,
												   double* reference_conformation,
												   int reference_conformation_id,
												   int number_of_conformations,
												   int number_of_atoms,
												   double* rmsd );

		virtual void calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoords(double* all_fit_coordinates,
																				   double* all_calc_coordinates,
																				   double*fit_reference,
																				   double*calc_reference,
																				   int reference_conformation_id,
																				   int number_of_conformations,
																				   int number_of_fit_atoms,
																				   int number_of_calc_atoms,
																				   double* rmsd);

		virtual void calcRMSDOfOneVsFollowingModifyingCoordinates( double* all_coordinates,
																		   double* reference_conformation,
																		   int reference_conformation_id,
																		   int number_of_conformations,
																		   int number_of_atoms,
																		   double* rmsd );

		virtual void calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoordsModifyingCoordinates(
				double* all_fit_coordinates,
				double* all_calc_coordinates,
				double*fit_reference,
				double*calc_reference,
				int reference_conformation_id,
				int number_of_conformations,
				int number_of_fit_atoms,
				int number_of_calc_atoms,
				double* rmsd  );


};
#endif
