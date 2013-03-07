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

QTRFITOmpKernel::QTRFITOmpKernel(int number_of_threads){
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

QTRFITOmpKernel::~QTRFITOmpKernel(){}

void QTRFITOmpKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(
		double* reference, int reference_conformation_number, double* rmsd,
		int numberOfConformations, int coordinatesPerConformation, int atomsPerConformation,
		double *allCoordinates) {

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < numberOfConformations;++i) {
		// Coordinates are copied so the working coordinates set is not modified
		double* conformation_coords = &(allCoordinates[i * coordinatesPerConformation]);
		double* conformation_coords_tmp = new double[coordinatesPerConformation];

		RMSDTools::copyArrays(conformation_coords_tmp, conformation_coords, coordinatesPerConformation);

		superpose(atomsPerConformation, conformation_coords_tmp, reference);
		rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(reference,
				conformation_coords_tmp, atomsPerConformation);
		delete[] conformation_coords_tmp;
	}
}

void QTRFITOmpKernel::oneVsFollowingFitEqualCalcWithConfRotation(
		double* reference, int reference_conformation_number, double* rmsd,
		int numberOfConformations, int coordinatesPerConformation, int atomsPerConformation,
		double *allCoordinates) {

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {
		// Real conformation coordinates are used, so they are modified
		double* conformation_coords = &(allCoordinates[i * coordinatesPerConformation]);

		superpose(atomsPerConformation, conformation_coords, reference);
		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					reference, conformation_coords, atomsPerConformation);
		}
	}
}

void QTRFITOmpKernel::oneVsFollowingFitDiffersCalcWithoutConfRotation(
		double* fitReference, double* calcReference, int reference_conformation_number,
		double* rmsd, int numberOfConformations,
		int coordinatesPerConformation, int atomsPerConformation, double *allCoordinates,
		int coordinatesPerRMSDConformation, int atomsPerRMSDConformation, double *allRMSDCoordinates){

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {

		int fitCoordsOffset = i * coordinatesPerConformation;
		double* fit_conformation_coords = &(allCoordinates[fitCoordsOffset]);
		double* fit_conformation_coords_tmp =	new double[coordinatesPerConformation];
		RMSDTools::copyArrays(fit_conformation_coords_tmp, fit_conformation_coords, coordinatesPerConformation);

		int calcCoordsOffset = i * coordinatesPerRMSDConformation;
		double* calc_conformation_coords =	&(allRMSDCoordinates[calcCoordsOffset]);
		double* calc_conformation_coords_tmp = 	new double[coordinatesPerRMSDConformation];
		RMSDTools::copyArrays(calc_conformation_coords_tmp, calc_conformation_coords, coordinatesPerRMSDConformation);

		superpose(atomsPerConformation, fit_conformation_coords_tmp,
				fitReference, atomsPerRMSDConformation,	calc_conformation_coords_tmp);

		rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
				calcReference, calc_conformation_coords_tmp, atomsPerRMSDConformation);

		delete[] fit_conformation_coords_tmp;
		delete[] calc_conformation_coords_tmp;
	}
}

void QTRFITOmpKernel::oneVsFollowingFitDiffersCalcWithConfRotation(
		double* fitReference, double* calcReference, int reference_conformation_number,
		double* rmsd, int numberOfConformations,
		int coordinatesPerConformation, int atomsPerConformation, double *allCoordinates,
		int coordinatesPerRMSDConformation, int atomsPerRMSDConformation, double *allRMSDCoordinates) {

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {

		double* fit_conformation_coords = &(allCoordinates[i* coordinatesPerConformation]);
		double* calc_conformation_coords = &(allRMSDCoordinates[i* coordinatesPerRMSDConformation]);

		superpose(atomsPerConformation, fit_conformation_coords,
				fitReference, atomsPerRMSDConformation,
				calc_conformation_coords);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					calcReference, calc_conformation_coords, atomsPerRMSDConformation);
		}
	}
}
