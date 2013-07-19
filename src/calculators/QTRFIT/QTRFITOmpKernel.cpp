/*
 * QTRFITOmpKernel.cpp
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#include "QTRFITOmpKernel.h"
#include "../RMSDTools.h"
#include <omp.h>
#include <cstddef>
#include "../RMSDCalculationData.h"

QTRFITOmpKernel::QTRFITOmpKernel(int number_of_threads){
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

QTRFITOmpKernel::~QTRFITOmpKernel(){}

void QTRFITOmpKernel::oneVsFollowingFitEqualCalcCoords(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < data->numberOfConformations;	++i) {
		// Real conformation coordinates are used, so they are modified
		double* conformation_coords = data->getFittingConformationAt(i);

		superpose(data->atomsPerFittingConformation,
						conformation_coords,
						reference);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL){
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					reference,
					conformation_coords,
					data->atomsPerFittingConformation);
		}
	}
}

void QTRFITOmpKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data) {

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < data->numberOfConformations;	++i) {

		double* fit_conformation_coords = data->getFittingConformationAt(i);
		double* calc_conformation_coords = data->getCalculationConformationAt(i);


		superpose(data->atomsPerFittingConformation,
						fit_conformation_coords,
						fitReference,
						data->atomsPerCalculationConformation,
						calc_conformation_coords);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}
