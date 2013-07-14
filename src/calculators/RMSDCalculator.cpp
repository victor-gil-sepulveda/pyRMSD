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
	this->atomsPerFittingConformation = atomsPerConformation;
	this->allFittingCoordinates = allCoordinates;
	this->coordinatesPerFittingConformation = atomsPerConformation*3;
	this->allCalculationCoordinates = NULL;
	this->coordinatesPerCalculationConformation = 0;
	this->atomsPerCalculationConformation = 0;
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
	this->atomsPerCalculationConformation = atomsPerRMSDConformation;
	this->coordinatesPerCalculationConformation = atomsPerRMSDConformation*3;
	this->allCalculationCoordinates = allRMSDCoordinates;

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
		if(this->allCalculationCoordinates == NULL){
			double* reference_conformation = &(allFittingCoordinates[reference_conformation_number*coordinatesPerFittingConformation]);

			this->_one_vs_following_fit_equals_calc_coords_rotating_coordinates(
					reference_conformation,
					reference_conformation_number,
					rmsd);
		}
		else{

			double* fit_reference_conformation = &(allFittingCoordinates[reference_conformation_number*coordinatesPerFittingConformation]);
			double* calc_reference_conformation = &(allCalculationCoordinates[reference_conformation_number*coordinatesPerCalculationConformation]);

			this->_one_vs_following_fit_differs_calc_coords_rotating_coordinates(
					fit_reference_conformation,
					calc_reference_conformation,
					reference_conformation_number,
					rmsd);
		}
	}
}

/*
 * This function creates the upper triangle of the rmsd matrix for a given coordinates set.
 * It doesn't modify coordinates at all.
 */
void RMSDCalculator::calculateRMSDCondensedMatrix(std::vector<double>& rmsd){

	int num_of_rmsds = numberOfConformations * (numberOfConformations-1.) /2.;
	double* rmsd_tmp =  new double[num_of_rmsds];

	// Coordinates are modified here
	double* fit_centers =new double[numberOfConformations*3];
	//double* calc_centers = NULL;
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, numberOfConformations, allFittingCoordinates, fit_centers);
	if(this->allCalculationCoordinates == NULL){
		//calc_centers = new double[numberOfConformations*3];
		//RMSDTools::centerAllAtOrigin(atomsPerCalculationConformation, numberOfConformations, allCalculationCoordinates, calc_centers);
		RMSDTools::applyTranslationsToAll(atomsPerCalculationConformation, numberOfConformations,
												allCalculationCoordinates, fit_centers, -1);
	}

	this->kernelFunctions->matrixInit(allFittingCoordinates,
										atomsPerFittingConformation*3,
										allCalculationCoordinates,
										atomsPerCalculationConformation*3,
										numberOfConformations);

	for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){
		int offset = conformation_number*(numberOfConformations-1)- (((conformation_number-1)*conformation_number)/2) ;
		double* fit_reference_conformation = &(allFittingCoordinates[conformation_number*coordinatesPerFittingConformation]);

		if(this->allCalculationCoordinates == NULL){
			this->kernelFunctions->matrixOneVsFollowingFitEqualCalcWithoutConfRotation(
																fit_reference_conformation,
																conformation_number,
																&(rmsd_tmp[offset]),
																numberOfConformations,
																coordinatesPerFittingConformation,
																atomsPerFittingConformation,
																allFittingCoordinates);
		}
		else{
			double* calc_reference_conformation = &(allCalculationCoordinates[conformation_number*coordinatesPerCalculationConformation]);
			this->kernelFunctions->matrixOneVsFollowingFitDiffersCalcWithoutConfRotation(
																fit_reference_conformation,
																calc_reference_conformation,
																conformation_number,
																&(rmsd_tmp[offset]),
																numberOfConformations,
																coordinatesPerFittingConformation,
																atomsPerFittingConformation,
																allFittingCoordinates,
																coordinatesPerCalculationConformation,
																atomsPerCalculationConformation,
																allCalculationCoordinates);
		}
	}

	this->kernelFunctions->matrixEnd(rmsd_tmp, num_of_rmsds, rmsd);
	delete [] rmsd_tmp;

	RMSDTools::applyTranslationsToAll(atomsPerFittingConformation, numberOfConformations, allFittingCoordinates, fit_centers);
	delete [] fit_centers;

	if(this->allCalculationCoordinates == NULL 	){
		RMSDTools::applyTranslationsToAll(atomsPerCalculationConformation, numberOfConformations, allCalculationCoordinates, fit_centers);
		//RMSDTools::applyTranslationsToAll(atomsPerCalculationConformation, numberOfConformations, allCalculationCoordinates, calc_centers);
		//delete [] calc_centers;
	}

}

void RMSDCalculator::_one_vs_following_fit_equals_calc_coords(
		double* reference,
		int reference_conformation_number,
		double *rmsd){
}

/*
 * This one modifies the input coordinates that will be superposed with the reference.
 *
 */
void RMSDCalculator::_one_vs_following_fit_equals_calc_coords_rotating_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, numberOfConformations, allFittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerFittingConformation,
			atomsPerFittingConformation,
			allFittingCoordinates);
}

void RMSDCalculator::_one_vs_following_fit_differs_calc_coords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double *rmsd){

}

void RMSDCalculator::_one_vs_following_fit_differs_calc_coords_rotating_coordinates(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double *rmsd){

	double* fitCenters = NULL;

	fitCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, numberOfConformations, allFittingCoordinates, fitCenters);
	// Move the calculation set towards the center but keeping its relative position with the fitting set
	RMSDTools::applyTranslationsToAll(atomsPerCalculationConformation, numberOfConformations, allCalculationCoordinates, fitCenters, -1);
	delete [] fitCenters;

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithConfRotation(
			fitReference,
			calcReference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerFittingConformation,
			atomsPerFittingConformation,
			allFittingCoordinates,
			coordinatesPerCalculationConformation,
			atomsPerCalculationConformation,
			allCalculationCoordinates);
}

void RMSDCalculator::iterativeSuperposition(double rmsd_diff_to_stop, double* iteration_rmsd){
	double MAX_ITERATIONS = 200;
	double* reference_coords = new double[coordinatesPerFittingConformation];
	double* mean_coords = new double[coordinatesPerFittingConformation];
	double rmsd_difference = 0.0;
	int current_iteration = 0;

	// In the first step, reference is the first conformation
	// reference = coordinates[0]
	RMSDTools::copyArrays(reference_coords, allFittingCoordinates, coordinatesPerFittingConformation);

	// Start iterative superposition
	do{
		rmsd_difference = iterativeSuperpositionStep(reference_coords, mean_coords);

		if (iteration_rmsd!=NULL){
			iteration_rmsd[current_iteration] = rmsd_difference;
		}

		current_iteration++;
	}
	while(rmsd_difference > rmsd_diff_to_stop and current_iteration < MAX_ITERATIONS);

	// One last superposition is performed, and the other rotation coordinates are moved here
	superposition_with_external_reference(reference_coords);

	// Some cleaning
	delete [] reference_coords;
	delete [] mean_coords;
}

double RMSDCalculator::iterativeSuperpositionStep(double* reference_coords, double* mean_coords ){
		double rmsd_difference = 0;

		superposition_with_external_reference(reference_coords);

		// Calculate new mean coords, which will be the next reference
		RMSDTools::calculateMeanCoordinates(mean_coords,
											allFittingCoordinates,
											numberOfConformations,
											atomsPerFittingConformation);

		// rmsd(reference,mean)
		rmsd_difference = RMSDTools::calcRMS(reference_coords,
												mean_coords,
												atomsPerFittingConformation);

		// reference = mean
		RMSDTools::copyArrays(reference_coords, mean_coords, coordinatesPerFittingConformation);

		return rmsd_difference;
}

void RMSDCalculator::superposition_with_external_reference(double* reference_coords){
	if (allCalculationCoordinates == NULL){
		superposition_with_external_reference_without_calc_coords(reference_coords);
	}
	else{
		superposition_with_external_reference_rotating_calc_coords(reference_coords);
	}
}


// Reference coordinates is already a copy, in a different memory space than allCoordinates
// In this case (used for iterative superposition) conformations are recentered over the reference conformation.
void RMSDCalculator::superposition_with_external_reference_without_calc_coords(double* reference){

	// Recentering coordinates each step gives us am easier way of comparing with prody,
	// and removes any translation artifact
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, 1, reference);
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, numberOfConformations, allFittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			-1,
			NULL,
			numberOfConformations,
			coordinatesPerFittingConformation,
			atomsPerFittingConformation,
			allFittingCoordinates);
}

void RMSDCalculator::superposition_with_external_reference_rotating_calc_coords(double* reference){

	// Center fitting coordinates
	double* fitCenters = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation, 1, reference);
	RMSDTools::centerAllAtOrigin(atomsPerFittingConformation,
									numberOfConformations,
									allFittingCoordinates,
									fitCenters);

	// Apply a relative translation to center also the coordinates we want to rotate
	// them with the fitting coordinates
	RMSDTools::applyTranslationsToAll(atomsPerCalculationConformation,
										numberOfConformations,
										allCalculationCoordinates,
										fitCenters,
										-1);
	delete [] fitCenters;

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithConfRotation(
			reference,
			NULL,
			-1,
			NULL,
			numberOfConformations,
			coordinatesPerFittingConformation,
			atomsPerFittingConformation,
			allFittingCoordinates,
			coordinatesPerCalculationConformation,
			atomsPerCalculationConformation,
			allCalculationCoordinates);

}
