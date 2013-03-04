/*
 * KernelFunctions.h
 *
 *  Created on: 01/03/2013
 *      Author: victor
 */

#ifndef QCP_KERNELFUNCTIONS_H_
#define QCP_KERNELFUNCTIONS_H_

#include <cstddef>
#include "../../KernelFunctions.h"

class QCPKernelFunctions:public KernelFunctions {
	public:
		QCPKernelFunctions(){}
		virtual ~QCPKernelFunctions(){};
		
		virtual double innerProduct(
				double* A,
				double* first_conformation_coords,
				double* second_conformation_coords,
				int number_of_atoms) = 0;

		virtual double calcRMSDForTwoConformationsWithTheobaldMethod(
				double *A,
				double E0,
				int number_of_atoms,
				double* rot_matrix = NULL) = 0;

		virtual double calcRMSDOfTwoConformations(
				double* first_conformation_coords,
				double* second_conformation_coords,
				int number_of_atoms,
				double* rot_matrix = NULL) = 0;

		virtual void calcRMSDOfOneVsFollowing(double* all_coordinates,
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

		virtual void calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoords(   double* all_fit_coordinates,
																					   double* all_calc_coordinates,
																					   double*fit_reference,
																					   double*calc_reference,
																					   int reference_conformation_id,
																					   int number_of_conformations,
																					   int number_of_fit_atoms,
																					   int number_of_calc_atoms,
																					   double* rmsd)=0;
		virtual void calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoordsModifyingCoordinates(
													double* all_fit_coordinates,
													double* all_calc_coordinates,
													double*fit_reference,
													double*calc_reference,
													int reference_conformation_id,
													int number_of_conformations,
													int number_of_fit_atoms,
													int number_of_calc_atoms,
													double* rmsd  ) = 0;

};

#endif /* QCP_KERNELFUNCTIONS_H_ */
