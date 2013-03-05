#include <cmath>
#include <iostream>
#include "kernel_functions_omp.h"
#include "../../RMSDTools.h"

using namespace std;

#define floating_point_type double

void ThRMSDSerialOmpKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(double* all_coordinates,
														   double* reference_conformation,
														   int reference_conformation_id,
														   int number_of_conformations,
														   int number_of_atoms,
														   double* rmsd){

	int coordinates_per_conformation = number_of_atoms * 3;

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = reference_conformation_id + 1;
					second_conformation_id < number_of_conformations; ++second_conformation_id){
		double* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);
		double rmsd_val = calcRMSDOfTwoConformations(	reference_conformation,
															  	second_conformation_coords,
																number_of_atoms);
		if(rmsd!= NULL){
			rmsd[second_conformation_id-(reference_conformation_id+1)] = rmsd_val;
		}
	}
}

void ThRMSDSerialOmpKernel::oneVsFollowingFitEqualCalcWithConfRotation( double* all_coordinates,
																					double* reference_conformation,
																					int reference_conformation_id,
																					int number_of_conformations,
																					int number_of_atoms,
																					double* rmsd ){
	int coordinates_per_conformation = number_of_atoms * 3;

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = reference_conformation_id + 1;
			second_conformation_id < number_of_conformations; ++second_conformation_id){

		double rot_matrix[9];

		double* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);

		RMSDTools::initializeTo(rot_matrix, 0.0, 9);

		double rmsd_val = calcRMSDOfTwoConformations(	reference_conformation,
													  	second_conformation_coords,
														number_of_atoms,
														rot_matrix);
		if(rmsd!=NULL){
			rmsd[second_conformation_id-(reference_conformation_id+1)] = rmsd_val;
		}

		RMSDTools::rotate3D(number_of_atoms, second_conformation_coords, rot_matrix);
	}
}
