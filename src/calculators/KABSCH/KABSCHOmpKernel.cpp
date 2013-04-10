/*
 * KABSCHOmpKernel.cpp
 *
 *  Created on: 07/03/2013
 *      Author: victor
 */

#include "KABSCHOmpKernel.h"
#include "../RMSDTools.h"
#include <cmath>
#include <iostream>
#include <omp.h>
using namespace std;

KABSCHOmpKernel::KABSCHOmpKernel(int number_of_threads) {
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

KABSCHOmpKernel::~KABSCHOmpKernel() {

}

void KABSCHOmpKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			int numberOfConformations,
			int coordinatesPerConformation,
			int atomsPerConformation,
			double *allCoordinates){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = reference_conformation_number + 1;
					second_conformation_id < numberOfConformations; ++second_conformation_id){
		double* second_conformation_coords = &(allCoordinates[second_conformation_id*coordinatesPerConformation]);
		double rmsd_val = this->calculate_rmsd(
				reference,
				second_conformation_coords,
				atomsPerConformation);

		rmsd[second_conformation_id-(reference_conformation_number+1)] = rmsd_val;
	}
}

void KABSCHOmpKernel::oneVsFollowingFitEqualCalcWithConfRotation(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			int numberOfConformations,
			int coordinatesPerConformation,
			int atomsPerConformation,
			double *allCoordinates){

		#pragma omp parallel for shared(rmsd)
		for (int second_conformation_id = reference_conformation_number + 1;
				second_conformation_id < numberOfConformations; ++second_conformation_id){

			double rot_matrix[3][3];
			double* second_conformation_coords = &(allCoordinates[second_conformation_id*coordinatesPerConformation]);

			RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

			double rmsd_val = this->calculate_rotation_rmsd(
					reference,
					second_conformation_coords,
					atomsPerConformation,
					rot_matrix);

			if(rmsd!=NULL){
				rmsd[second_conformation_id-(reference_conformation_number+1)] = rmsd_val;
			}

			RMSDTools::rotate3D(atomsPerConformation, second_conformation_coords, rot_matrix);
		}
}

void KABSCHOmpKernel::oneVsFollowingFitDiffersCalcWithoutConfRotation(
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
		double *allRMSDCoordinates){

	int coordinates_per_fit_conformation = atomsPerConformation * 3;
	int coordinates_per_calc_conformation = atomsPerRMSDConformation * 3;

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = reference_conformation_number + 1;
			second_conformation_id < numberOfConformations; ++second_conformation_id){

		double rot_matrix[3][3];
		double* fit_conformation_coords = &(allCoordinates[second_conformation_id*coordinates_per_fit_conformation]);
		RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

		this->calculate_rotation_rmsd(
				fitReference,
				fit_conformation_coords,
				atomsPerConformation,
				rot_matrix);

		double* calc_conformation_coords =  &(allRMSDCoordinates[second_conformation_id*coordinates_per_calc_conformation]);
		double* calc_conformation_coords_copy = new double[coordinates_per_calc_conformation];
		RMSDTools::copyArrays(calc_conformation_coords_copy,calc_conformation_coords,coordinates_per_calc_conformation);
		RMSDTools::rotate3D(atomsPerRMSDConformation, calc_conformation_coords_copy, rot_matrix);

		rmsd[second_conformation_id-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
				calc_conformation_coords_copy, atomsPerRMSDConformation);

		delete [] calc_conformation_coords_copy;
	}
}

void KABSCHOmpKernel::oneVsFollowingFitDiffersCalcWithConfRotation(
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
		double *allRMSDCoordinates){

	int coordinates_per_fit_conformation = atomsPerConformation * 3;
	int coordinates_per_calc_conformation = atomsPerRMSDConformation * 3;

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = reference_conformation_number + 1;
			second_conformation_id < numberOfConformations; ++second_conformation_id){

		double rot_matrix[3][3];
		double* fit_conformation_coords = &(allCoordinates[second_conformation_id*coordinates_per_fit_conformation]);
		RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

		this->calculate_rotation_rmsd(
						fitReference,
						fit_conformation_coords,
						atomsPerConformation,
						rot_matrix);


		double* calc_conformation_coords =  &(allRMSDCoordinates[second_conformation_id*coordinates_per_calc_conformation]);

		RMSDTools::rotate3D(atomsPerConformation, fit_conformation_coords, rot_matrix);
		RMSDTools::rotate3D(atomsPerRMSDConformation, calc_conformation_coords, rot_matrix);

		if(rmsd!=NULL){
			rmsd[second_conformation_id-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
					calc_conformation_coords, atomsPerRMSDConformation);
		}
	}
}

