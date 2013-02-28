/*
 * RMSDomp.cpp
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#include "RMSDomp.h"
#include "../RMSDTools.h"
#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
using namespace std;


RMSDomp::RMSDomp(int numberOfConformations, int atomsPerConformation, double *allCoordinates):
		RMSD(numberOfConformations, atomsPerConformation,allCoordinates){
}

RMSDomp::~RMSDomp(){}

void RMSDomp::_one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	// Indeed only numberOfConformations-reference_conformation_number+1 conformations have to be centered...
	// TODO: change this to gain performance
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number+1; i < numberOfConformations; ++i){

		// Coordinates are copied so the working coordinates set is not modified
		double* conformation_coords = &(allCoordinates[i*coordinatesPerConformation]);
		double* conformation_coords_tmp = new double[coordinatesPerConformation];
		copy(conformation_coords,conformation_coords+coordinatesPerConformation,conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, conformation_coords_tmp, reference);

		rmsd[i-(reference_conformation_number+1)] = RMSDTools::calcRMS(reference, conformation_coords_tmp, atomsPerConformation);

		delete [] conformation_coords_tmp;
	}

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, centers);
	delete [] centers;
}

void RMSDomp::_one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
	delete [] centers;

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number+1; i < numberOfConformations; ++i){
		// Real conformation coordinates are used, so they are modified
		double* conformation_coords = &(allCoordinates[i*coordinatesPerConformation]);

		RMSDTools::superpose(atomsPerConformation, conformation_coords, reference);
		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL){
			rmsd[i-(reference_conformation_number+1)] = RMSDTools::calcRMS(reference, conformation_coords, atomsPerConformation);
		}
	}
	// Coordinates are left centered
}

void RMSDomp::_one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number+1; i < numberOfConformations; ++i){
		int fitCoordsOffset = i*coordinatesPerConformation;
		double* fit_conformation_coords = &(allCoordinates[fitCoordsOffset]);
		double* fit_conformation_coords_tmp = new double[coordinatesPerConformation];
		copy(fit_conformation_coords,fit_conformation_coords+coordinatesPerConformation,fit_conformation_coords_tmp);

		int calcCoordsOffset = i*coordinatesPerRMSDConformation;
		double* calc_conformation_coords = &(allRMSDCoordinates[calcCoordsOffset]);
		double* calc_conformation_coords_tmp = new double[coordinatesPerRMSDConformation];
		copy(calc_conformation_coords, calc_conformation_coords+coordinatesPerRMSDConformation, calc_conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, 	   fit_conformation_coords_tmp,  fitReference,
							 atomsPerRMSDConformation, calc_conformation_coords_tmp, calcReference);

		rmsd[i-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference, calc_conformation_coords_tmp, atomsPerRMSDConformation);

		delete [] fit_conformation_coords_tmp;
		delete [] calc_conformation_coords_tmp;
	}

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, fitCenters);
	RMSDTools::applyTranslationsToAll(this->atomsPerRMSDConformation, this->numberOfConformations, this->allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;
}

void RMSDomp::_one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;

	#pragma omp parallel for shared(rmsd)
	for (int i = reference_conformation_number+1; i < numberOfConformations; ++i){

		double* fit_conformation_coords = &(allCoordinates[i*coordinatesPerConformation]);
		double* calc_conformation_coords = &(allRMSDCoordinates[i*coordinatesPerRMSDConformation]);

		RMSDTools::superpose(atomsPerConformation, 	   fit_conformation_coords,  fitReference,
							 atomsPerRMSDConformation, calc_conformation_coords, calcReference);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL){
			rmsd[i-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference, calc_conformation_coords, atomsPerRMSDConformation);
		}

	}

}

