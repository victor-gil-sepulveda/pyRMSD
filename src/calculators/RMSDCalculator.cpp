#include "RMSDCalculator.h"
#include "RMSDTools.h"
#include <iostream>
#include <cstdlib>
#include <vector>
using namespace std;

RMSDCalculator::RMSDCalculator(	int numberOfConformations,
								int atomsPerConformation,
								double* allCoordinates,
								KernelFunctions* kernelFunctions){
	this->numberOfConformations = numberOfConformations;
	this->atomsPerConformation = atomsPerConformation;
	this->allCoordinates = allCoordinates;
	this->coordinatesPerConformation = atomsPerConformation*3;
	this->allRMSDCoordinates = NULL;
	this->coordinatesPerRMSDConformation = 0;
	this->atomsPerRMSDConformation = 0;
	this->rotateFittingCoordinates = false;
	this->kernelFunctions = kernelFunctions;
}

RMSDCalculator::~RMSDCalculator(){
	if(this->kernelFunctions != NULL){
		delete this->kernelFunctions;
	}
}

/*
 * Sets a different set of coordinates for RMSD calculation and fit.
 */
void RMSDCalculator::setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates){
	this->atomsPerRMSDConformation = atomsPerRMSDConformation;
	this->coordinatesPerRMSDConformation = atomsPerRMSDConformation*3;
	this->allRMSDCoordinates = allRMSDCoordinates;
	this->kernelFunctions->changeCalculationCoords(
			allRMSDCoordinates,
			atomsPerRMSDConformation,
			numberOfConformations);
}

void RMSDCalculator::oneVsFollowing(int reference_conformation_number, double* rmsd){

	if(reference_conformation_number >= numberOfConformations){
		cout<<"Error, this conformation doesn't exist ("<<reference_conformation_number<<")"<<endl;
	}
	else{
		if(this->allRMSDCoordinates == NULL){
			double* reference_conformation = &(allCoordinates[reference_conformation_number*coordinatesPerConformation]);

			if(!this->rotateFittingCoordinates){
				this->_one_vs_following_fit_equals_calc_coords(
						reference_conformation,
						reference_conformation_number,
						rmsd);
			}
			else{
				this->_one_vs_following_fit_equals_calc_coords_rotating_coordinates(
						reference_conformation,
						reference_conformation_number,
						rmsd);
			}
		}
		else{
			double* fit_reference_conformation = &(allCoordinates[reference_conformation_number*coordinatesPerConformation]);
			double* calc_reference_conformation = &(allRMSDCoordinates[reference_conformation_number*coordinatesPerRMSDConformation]);

			if(!this->rotateFittingCoordinates){
				this->_one_vs_following_fit_differs_calc_coords(
						fit_reference_conformation,
						calc_reference_conformation,
						reference_conformation_number,
						rmsd);
			}
			else{
				this->_one_vs_following_fit_differs_calc_coords_rotating_coordinates(
						fit_reference_conformation,
						calc_reference_conformation,
						reference_conformation_number,
						rmsd);
			}
		}
	}
}

void RMSDCalculator::calculateRMSDCondensedMatrix(std::vector<double>& rmsd){

	int num_of_rmsds = numberOfConformations * (numberOfConformations-1.) /2.;
	double* rmsd_tmp =  new double[num_of_rmsds];

	// Coordinates are modified here
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates);

	for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){
		int offset = conformation_number*(numberOfConformations-1)- (((conformation_number-1)*conformation_number)/2) ;

		if(this->allRMSDCoordinates == NULL){
			double* reference_conformation = &(allCoordinates[conformation_number*coordinatesPerConformation]);

			this->kernelFunctions->oneVsFollowingFitEqualCalcWithoutConfRotation(
																reference_conformation,
																conformation_number,
																&(rmsd_tmp[offset]),
																numberOfConformations,
																coordinatesPerConformation,
																atomsPerConformation,
																allCoordinates);
		}
		else{
			double* fit_reference_conformation = &(allCoordinates[conformation_number*coordinatesPerConformation]);
			double* calc_reference_conformation = &(allRMSDCoordinates[conformation_number*coordinatesPerRMSDConformation]);

			this->kernelFunctions->oneVsFollowingFitDiffersCalcWithoutConfRotation(
																fit_reference_conformation,
																calc_reference_conformation,
																conformation_number,
																&(rmsd_tmp[offset]),
																numberOfConformations,
																coordinatesPerConformation,
																atomsPerConformation,
																allCoordinates,
																coordinatesPerRMSDConformation,
																atomsPerRMSDConformation,
																allRMSDCoordinates);
		}
	}

	for (int i = 0; i < num_of_rmsds; ++i){
		rmsd.push_back(rmsd_tmp[i]);
	}

	delete [] rmsd_tmp;
}

void RMSDCalculator::_one_vs_following_fit_equals_calc_coords(
		double* reference,
		int reference_conformation_number,
		double *rmsd){


	double* centers = NULL;

	centers = new double[(numberOfConformations-reference_conformation_number)*3]; // A bit of memory is wasted
	RMSDTools::centerAllAtOrigin(
			atomsPerConformation,
			(numberOfConformations-reference_conformation_number),
			reference,
			centers);

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithoutConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates);

	// Move then again to their places
	RMSDTools::applyTranslationsToAll(
			this->atomsPerConformation,
			(numberOfConformations-reference_conformation_number),
			reference,
			centers);
	delete [] centers;
}

void RMSDCalculator::_one_vs_following_fit_equals_calc_coords_rotating_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
	delete [] centers;

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates);
}

void RMSDCalculator::_one_vs_following_fit_differs_calc_coords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double *rmsd){

	double* fitCenters = NULL;
	double* calcCenters = NULL;

	fitCenters = new double[numberOfConformations*3];
	calcCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithoutConfRotation(
			fitReference,
			calcReference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates,
			coordinatesPerRMSDConformation,
			atomsPerRMSDConformation,
			allRMSDCoordinates);

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, fitCenters);
	RMSDTools::applyTranslationsToAll(this->atomsPerRMSDConformation, this->numberOfConformations, this->allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;
}

void RMSDCalculator::_one_vs_following_fit_differs_calc_coords_rotating_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates);

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithConfRotation(
			fitReference,
			calcReference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates,
			coordinatesPerRMSDConformation,
			atomsPerRMSDConformation,
			allRMSDCoordinates);
}

void RMSDCalculator::iterativeSuperposition(double rmsd_diff_to_stop){
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
		if(this->allRMSDCoordinates == NULL){
			superposition_with_external_reference_and_fit_equals_calc(reference_coords, NULL);
		}
		else{
			// As we are only interested in coordinates change, we can set the calc_reference to NULL
			// (we won't calculate the RMSD)
			superposition_with_external_reference_and_fit_differs_calc(reference_coords);
		}
		// Calculate new mean coords, which will be the next reference
		RMSDTools::calculateMeanCoordinates(mean_coords, allCoordinates,
											numberOfConformations, atomsPerConformation);

		// rmsd(reference,mean)
		rmsd_difference = RMSDTools::calcRMS(reference_coords, mean_coords, atomsPerConformation);

		// reference = mean
		RMSDTools::copyArrays(reference_coords, mean_coords, coordinatesPerConformation);

		//cout<<"step: "<< current_iteration<<" rmsd diff: "<<rmsd_difference<<endl;

		current_iteration++;
	}
	while(rmsd_difference > rmsd_diff_to_stop and	current_iteration < MAX_ITERATIONS);

	delete [] reference_coords;
	delete [] mean_coords;
}

// Reference coordinates is already a copy, in a different memory space than allCoordinates
// In this case (used for iterative superposition) conformations are recentered over the reference conformation.
void RMSDCalculator::superposition_with_external_reference_and_fit_equals_calc(double* reference, double* rmsds){
	double reference_center[3];

	RMSDTools::centerAllAtOrigin(atomsPerConformation, 1, reference, reference_center);

	_one_vs_following_fit_equals_calc_coords_rotating_coordinates(reference, -1, rmsds);

	RMSDTools::applyTranslationToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, reference_center);
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, 1, reference, reference_center);
}

void RMSDCalculator::superposition_with_external_reference_and_fit_differs_calc(double* fit_reference){
	double fit_reference_center[3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, 1, fit_reference, fit_reference_center);

	_one_vs_following_fit_differs_calc_coords_rotating_coordinates(fit_reference, NULL, -1, NULL);

	RMSDTools::applyTranslationToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, fit_reference_center);
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, 1, fit_reference, fit_reference_center);
}
