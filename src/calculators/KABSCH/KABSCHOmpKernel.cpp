/*
 * KABSCHOmpKernel.cpp
 *
 *  Created on: 07/03/2013
 *      Author: victor
 */

#include "KABSCHOmpKernel.h"
#include "../RMSDTools.h"
#include <cmath>
#include <omp.h>
#include "../RMSDCalculationData.h"

KABSCHOmpKernel::KABSCHOmpKernel(int number_of_threads) {
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

KABSCHOmpKernel::~KABSCHOmpKernel() {

}

void KABSCHOmpKernel::oneVsFollowingFitEqualCalcCoords(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			RMSDCalculationData* data){

		#pragma omp parallel for shared(rmsd)
		for (int second_conformation_index = reference_conformation_number + 1;
				second_conformation_index < data->numberOfConformations; ++second_conformation_index){

			double rot_matrix[3][3];

			double* second_conformation_coords = data->getFittingConformationAt(second_conformation_index);

			RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

			double rmsd_val = this->calculate_rotation_rmsd(
					reference,
					second_conformation_coords,
					data->atomsPerFittingConformation,
					rot_matrix);

			if(rmsd!=NULL){
				rmsd[second_conformation_index-(reference_conformation_number+1)] = rmsd_val;
			}

			RMSDTools::rotate3D(data->atomsPerFittingConformation,
								second_conformation_coords,
								rot_matrix);
		}
}

void KABSCHOmpKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double rot_matrix[3][3];

		double* fit_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		double* calc_conformation_coords =  data->getCalculationConformationAt(second_conformation_index);
		RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

		this->calculate_rotation_rmsd(
						fitReference,
						fit_conformation_coords,
						data->atomsPerFittingConformation,
						rot_matrix);

		RMSDTools::rotate3D(data->atomsPerFittingConformation,
				fit_conformation_coords,
				rot_matrix);

		RMSDTools::rotate3D(data->atomsPerCalculationConformation,
				calc_conformation_coords,
				rot_matrix);

		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}

