/*
 * KernelFunctions.h
 *
 *  Created on: 01/03/2013
 *      Author: victor
 */

#ifndef KERNELFUNCTIONS_H_
#define KERNELFUNCTIONS_H_
#include <cstddef>

class KernelFunctions {
	public:
		KernelFunctions(){}
		virtual ~KernelFunctions(){};
		
		virtual double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms) = 0;

		virtual double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms, double* rot_matrix = NULL) = 0;

		virtual double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords,
				int number_of_atoms, double* rot_matrix = NULL) = 0;

		virtual void calcRMSDOfOneVsFollowing( double* all_coordinates,
												   double* reference_conformation,
												   int reference_conformation_id,
												   int number_of_conformations,
												   int number_of_atoms,
												   double* rmsd ) = 0;

		virtual void calcRMSDOfOneVsFollowingModifyingCoordinates( double* all_coordinates,
																		   double* reference_conformation,
																		   int reference_conformation_id,
																		   int number_of_conformations,
																		   int number_of_atoms,
																		   double* rmsd ) = 0;
};

#endif /* KERNELFUNCTIONS_H_ */
