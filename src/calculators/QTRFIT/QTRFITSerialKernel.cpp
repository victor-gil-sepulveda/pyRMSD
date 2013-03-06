/*
 * QTRFITSerialKernel.cpp
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#include "QTRFITSerialKernel.h"
#include <cstddef>
#include <algorithm>
#include "../RMSDTools.h"

using std::copy;

QTRFITSerialKernel::QTRFITSerialKernel() {}

QTRFITSerialKernel::~QTRFITSerialKernel() {}

void QTRFITSerialKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(
		double* reference, int reference_conformation_number, double* rmsd,
		int numberOfConformations, int coordinatesPerConformation, int atomsPerConformation,
		double *allCoordinates) {

	for (int i = reference_conformation_number + 1; i < numberOfConformations;++i) {
		// Coordinates are copied so the working coordinates set is not modified
		double* conformation_coords = &(allCoordinates[i * coordinatesPerConformation]);
		double* conformation_coords_tmp = new double[coordinatesPerConformation];

		copy(conformation_coords, conformation_coords + coordinatesPerConformation, conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, conformation_coords_tmp, reference);
		rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(reference,
				conformation_coords_tmp, atomsPerConformation);
		delete[] conformation_coords_tmp;
	}
}

void QTRFITSerialKernel::oneVsFollowingFitEqualCalcWithConfRotation(
		double* reference, int reference_conformation_number, double* rmsd,
		int numberOfConformations, int coordinatesPerConformation, int atomsPerConformation,
		double *allCoordinates) {

	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {
		// Real conformation coordinates are used, so they are modified
		double* conformation_coords = &(allCoordinates[i * coordinatesPerConformation]);

		RMSDTools::superpose(atomsPerConformation, conformation_coords, reference);
		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					reference, conformation_coords, atomsPerConformation);
		}
	}
}

void QTRFITSerialKernel::oneVsFollowingFitDiffersCalcWithoutConfRotation(
		double* fitReference, double* calcReference, int reference_conformation_number,
		double* rmsd, int numberOfConformations,
		int coordinatesPerConformation, int atomsPerConformation, double *allCoordinates,
		int coordinatesPerRMSDConformation, int atomsPerRMSDConformation, double *allRMSDCoordinates){

	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {

		int fitCoordsOffset = i * coordinatesPerConformation;
		double* fit_conformation_coords = &(allCoordinates[fitCoordsOffset]);
		double* fit_conformation_coords_tmp =	new double[coordinatesPerConformation];
		copy(fit_conformation_coords, fit_conformation_coords + coordinatesPerConformation,
				fit_conformation_coords_tmp);

		int calcCoordsOffset = i * coordinatesPerRMSDConformation;
		double* calc_conformation_coords =	&(allRMSDCoordinates[calcCoordsOffset]);
		double* calc_conformation_coords_tmp = 	new double[coordinatesPerRMSDConformation];
		copy(calc_conformation_coords, calc_conformation_coords + coordinatesPerRMSDConformation,
				calc_conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, fit_conformation_coords_tmp,
				fitReference, atomsPerRMSDConformation,	calc_conformation_coords_tmp);

		rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
				calcReference, calc_conformation_coords_tmp, atomsPerRMSDConformation);

		delete[] fit_conformation_coords_tmp;
		delete[] calc_conformation_coords_tmp;
	}
}

void QTRFITSerialKernel::oneVsAllFitDiffersCalcWithConfRotation(
		double* fitReference, double* calcReference, int reference_conformation_number,
		double* rmsd, int numberOfConformations,
		int coordinatesPerConformation, int atomsPerConformation, double *allCoordinates,
		int coordinatesPerRMSDConformation, int atomsPerRMSDConformation, double *allRMSDCoordinates) {

	for (int i = reference_conformation_number + 1; i < numberOfConformations;	++i) {

		double* fit_conformation_coords = &(allCoordinates[i* coordinatesPerConformation]);
		double* calc_conformation_coords = &(allRMSDCoordinates[i* coordinatesPerRMSDConformation]);

		RMSDTools::superpose(atomsPerConformation, fit_conformation_coords,
				fitReference, atomsPerRMSDConformation,
				calc_conformation_coords);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					calcReference, calc_conformation_coords, atomsPerRMSDConformation);
		}
	}
}
