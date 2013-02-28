/*
 * RMSDomp.cpp
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#include "RMSDomp.h"
#include "RMSDTools.h"
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

RMSDomp::~RMSDomp(){
}

void RMSDomp::calculateRMSDCondensedMatrix(std::vector<double> & rmsd){
	int num_of_rmsds = numberOfConformations * (numberOfConformations-1.) /2.;
	double* rmsd_tmp =  new double[num_of_rmsds];

	for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){
		int offset = conformation_number*(numberOfConformations-1)- (((conformation_number-1)*conformation_number)/2) ;
		oneVsFollowing(conformation_number,&(rmsd_tmp[offset]));
	}

	for (int i = 0; i < num_of_rmsds; ++i){
		rmsd.push_back(rmsd_tmp[i]);
	}

	delete [] rmsd_tmp;
}

void RMSDomp::oneVsFollowing(int reference_conformation_number, double *rmsd){
	if(reference_conformation_number >= numberOfConformations){
		cout<<"Error, this conformation doesn't exist ("<<reference_conformation_number<<")"<<endl;
	}
	else{

		if(this->allRMSDCoordinates == NULL){
			double* reference_conformation = &(allCoordinates[reference_conformation_number*coordinatesPerConformation]);

			if(!this->modifyFittingCoordinates){
				this->_one_vs_following_fit_equals_calc_coords(reference_conformation, reference_conformation_number, rmsd);
			}
			else{
				this->_one_vs_following_fit_equals_calc_coords_changing_coordinates(reference_conformation, reference_conformation_number, rmsd);
			}
		}
		else{
			double* fit_reference_conformation = &(allCoordinates[reference_conformation_number*coordinatesPerConformation]);
			double* calc_reference_conformation = &(allRMSDCoordinates[reference_conformation_number*coordinatesPerRMSDConformation]);

			if(!this->modifyFittingCoordinates){
				this->_one_vs_following_fit_differs_calc_coords(fit_reference_conformation, calc_reference_conformation, reference_conformation_number, rmsd);
			}
			else{
				this->_one_vs_following_fit_differs_calc_coords_changing_coordinates(fit_reference_conformation, calc_reference_conformation, reference_conformation_number, rmsd);
			}
		}
	}
}

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

// Reference coords is already a copy, in a different memory space than allCoordinates
// In this case (used for iterative superposition) conformations are recentered over the reference conformation.
void RMSDomp::superposition_with_external_reference(double* reference, double* rmsds){
	double reference_center[3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, 1, reference, reference_center);

	_one_vs_following_fit_equals_calc_coords_changing_coordinates(reference, -1, rmsds);

	RMSDTools::applyTranslationToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, reference_center);
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, 1, reference, reference_center);
}

void RMSDomp::iterativeSuperposition(double rmsd_diff_to_stop){
	double MAX_ITERATIONS = 200;
	double* reference_coords = new double[coordinatesPerConformation];
	double* mean_coords = new double[coordinatesPerConformation];
	double rmsd_difference = 0.0;
	int current_iteration = 0;

	// In the first step, reference is the first conformation
	// reference = coordinates[0]
	RMSDTools::copyArrays(reference_coords, allCoordinates, coordinatesPerConformation);
	do{
		// Superpose all conformations with the reference conformation
		superposition_with_external_reference(reference_coords, NULL);

		// Calculate new mean coords, which will be the next reference
		RMSDTools::calculateMeanCoordinates(mean_coords, allCoordinates,
											numberOfConformations, atomsPerConformation);

		// rmsd(reference,mean)
		rmsd_difference = RMSDTools::calcRMS(reference_coords, mean_coords, atomsPerConformation);

		// reference = mean
		RMSDTools::copyArrays(reference_coords, mean_coords, coordinatesPerConformation);

		cout<<"step: "<< current_iteration<<" rmsd diff: "<<rmsd_difference<<endl;

		current_iteration++;
	}
	while(rmsd_difference > rmsd_diff_to_stop and
			current_iteration < MAX_ITERATIONS);
}
