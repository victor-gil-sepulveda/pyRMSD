#include "RMSDCalculator.h"
#include "RMSDTools.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include "RMSDCalculationData.h"
using namespace std;


/**
 * Class constructor. It receives a data package for the RMSD calculation and the implementation
 * of the RMSD calculation functions into the KernelFucntions object. This two objects are from now
 * on responsibility of this RMSDCalculator instance and cannot be shared
 *
 *  \param rmsdData Data package with coordinates and lengths.
 *  \param kernelFunctions Kernel object with RMSD calculation functions for a concrete method and
 *  parallelization scheme.
 *
 */
RMSDCalculator::RMSDCalculator(RMSDCalculationData* rmsdData, KernelFunctions* kernelFunctions){
	this->rmsdData = rmsdData;
	this->kernelFunctions = kernelFunctions;
	if(this->rmsdData->hasCalculationCoordinatesSet()){
		this->kernelFunctions->changeCalculationCoords(
				this->rmsdData->calculationCoordinates,
				this->rmsdData->atomsPerCalculationConformation,
				this->rmsdData->numberOfConformations);
	}
}

/**
 * Class destructor. It is responsible of destroying the data and kernel objects.
 */
RMSDCalculator::~RMSDCalculator(){
	delete this->rmsdData;
	delete this->kernelFunctions;
}

/**
 * 	Calculates the RMSD of one conformation VS the rest of conformations with greater index. All
 * 	input coordinates that are used get centered and superposed to the reference. This means that
 * 	if the index is > 0, conformations in 0:i are not modified.
 *
 *  \param reference_conformation_index One index in the range 0:N-1 where N is the total number
 *  of conformations.
 *  \param rmsd A pointer to an array with enough room to hold the calculated rmsds.
 */
void RMSDCalculator::oneVsFollowing(int reference_conformation_index, double* rmsd){

	if(reference_conformation_index >= this->rmsdData->numberOfConformations){
		cout<<"Error, this conformation doesn't exist ("<<reference_conformation_index<<")"<<endl;
	}
	else{
		if(!this->rmsdData->hasCalculationCoordinatesSet()){
			double* reference_conformation = this->rmsdData->getFittingConformationAt(reference_conformation_index);
			this->_one_vs_following_fit_equals_calc_coords_rotating_coordinates(
					reference_conformation,
					reference_conformation_index,
					rmsd);
		}
		else{

			double* fit_reference_conformation = this->rmsdData->getFittingConformationAt(reference_conformation_index);
			double* calc_reference_conformation = this->rmsdData->getCalculationConformationAt(reference_conformation_index);

			this->_one_vs_following_fit_differs_calc_coords_rotating_coordinates(
					fit_reference_conformation,
					calc_reference_conformation,
					reference_conformation_index,
					rmsd);
		}
	}
}

/**
 * Calculates the upper triangle of the RMSD matrix for a given coordinates set.
 * Coordinates get modified (iteratively centered and superposed).
 *
 * \param rmsd Linear vector that will hold the rmsd values of the matrix. The Nth-1 elements will
 * form the first row, the next N-2 the second row and so on.
 */
void RMSDCalculator::calculateRMSDCondensedMatrix(std::vector<double>& rmsd){

	int num_of_rmsds = this->rmsdData->numberOfConformations *
						(this->rmsdData->numberOfConformations-1.) /2.;
	rmsd.clear();
	rmsd.resize(num_of_rmsds);

	if(! this->rmsdData->hasCalculationCoordinatesSet()){
		calculateRMSDCondensedMatrixWithFittingCoordinates(rmsd);
	}
	else{
		calculateRMSDCondensedMatrixWithFittingAndCalculationCoordinates(rmsd);
	}
}

/**
 * Does the actual calculation of the matrix, doing the RMSD calculation using the
 * same coordinates that were used for fitting.
 *
 * \param rmsd Linear vector that will hold the rmsd values of the matrix. The Nth-1 elements will
 * form the first row, the next N-2 the second row and so on.
 */
void RMSDCalculator::calculateRMSDCondensedMatrixWithFittingCoordinates(vector<double>& rmsd){

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
									this->rmsdData->numberOfConformations,
									this->rmsdData->fittingCoordinates);

	this->kernelFunctions->matrixInit(this->rmsdData->fittingCoordinates,
										this->rmsdData->fittingConformationLength,
										this->rmsdData->calculationCoordinates,
										this->rmsdData->calculationConformationLength,
										this->rmsdData->numberOfConformations);

	for(int reference_index = 0; reference_index<this->rmsdData->numberOfConformations; ++reference_index){
			int offset = reference_index*(this->rmsdData->numberOfConformations-1)- (((reference_index-1)*reference_index)/2) ;
			double* reference_conformation = this->rmsdData->getFittingConformationAt(reference_index);
			this->kernelFunctions->matrixOneVsFollowingFitEqualCalcWithoutConfRotation(
					reference_conformation,
					reference_index,
					&(rmsd[offset]),
					this->rmsdData->numberOfConformations,
					this->rmsdData->fittingConformationLength,
					this->rmsdData->atomsPerFittingConformation,
					this->rmsdData->fittingCoordinates);
	}

}
/**
 * Exactly the same than 'calculateRMSDCondensedMatrixWithFittingCoordinates' but this time uses
 * one coordinate set for fitting and other to calculate the RMSD.
 *
 * \param rmsd Linear vector that will hold the rmsd values of the matrix. The Nth-1 elements will
 * form the first row, the next N-2 the second row and so on.
 */
void RMSDCalculator::calculateRMSDCondensedMatrixWithFittingAndCalculationCoordinates(vector<double>& rmsd){
	double* fit_centers =new double[this->rmsdData->numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
									this->rmsdData->numberOfConformations,
									this->rmsdData->fittingCoordinates,
									fit_centers);

	RMSDTools::applyTranslationsToAll(this->rmsdData->atomsPerCalculationConformation,
											this->rmsdData->numberOfConformations,
											this->rmsdData->calculationCoordinates,
											fit_centers,
											-1);

	delete [] fit_centers;

	this->kernelFunctions->matrixInit(this->rmsdData->fittingCoordinates,
										this->rmsdData->fittingConformationLength,
										this->rmsdData->calculationCoordinates,
										this->rmsdData->calculationConformationLength,
										this->rmsdData->numberOfConformations);

	for(int reference_index = 0; reference_index<this->rmsdData->numberOfConformations; ++reference_index){
		int offset = reference_index*(this->rmsdData->numberOfConformations-1)- (((reference_index-1)*reference_index)/2) ;
		double* fit_reference_conformation = this->rmsdData->getFittingConformationAt(reference_index);
		double* calc_reference_conformation = this->rmsdData->getCalculationConformationAt(reference_index);
		this->kernelFunctions->matrixOneVsFollowingFitDiffersCalcWithoutConfRotation(
															fit_reference_conformation,
															calc_reference_conformation,
															reference_index,
															&(rmsd[offset]),
															this->rmsdData->numberOfConformations,
															this->rmsdData->fittingConformationLength,
															this->rmsdData->atomsPerFittingConformation,
															this->rmsdData->fittingCoordinates,
															this->rmsdData->calculationConformationLength,
															this->rmsdData->atomsPerCalculationConformation,
															this->rmsdData->calculationCoordinates);
	}

}
/*
 *
 */
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
	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingConformationLength,
			this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->fittingCoordinates);
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

	fitCenters = new double[this->rmsdData->numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates,
			fitCenters);
	// Move the calculation set towards the center but keeping its relative position with the fitting set
	RMSDTools::applyTranslationsToAll(this->rmsdData->atomsPerCalculationConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->calculationCoordinates,
			fitCenters,
			-1);
	delete [] fitCenters;

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithConfRotation(
			fitReference,
			calcReference,
			reference_conformation_number,
			rmsd,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingConformationLength,
			this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->fittingCoordinates,
			this->rmsdData->calculationConformationLength,
			this->rmsdData->atomsPerCalculationConformation,
			this->rmsdData->calculationCoordinates);
}

void RMSDCalculator::iterativeSuperposition(double rmsd_diff_to_stop, double* iteration_rmsd){
	double MAX_ITERATIONS = 200;
	double* reference_coords = new double[this->rmsdData->fittingConformationLength];
	double* mean_coords = new double[this->rmsdData->fittingConformationLength];
	double rmsd_difference = 0.0;
	int current_iteration = 0;

	// In the first step, reference is the first conformation
	// reference = coordinates[0]
	RMSDTools::copyArrays(reference_coords,
			this->rmsdData->fittingCoordinates,
			this->rmsdData->fittingConformationLength);

	// Start iterative superposition
	do{
		rmsd_difference = iterative_superposition_step(reference_coords, mean_coords);

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

double RMSDCalculator::iterative_superposition_step(double* reference_coords, double* mean_coords ){
		double rmsd_difference = 0;

		superposition_with_external_reference(reference_coords);

		// Calculate new mean coords, which will be the next reference
		RMSDTools::calculateMeanCoordinates(mean_coords,
				this->rmsdData->fittingCoordinates,
				this->rmsdData->numberOfConformations,
				this->rmsdData->atomsPerFittingConformation);

		// rmsd(reference,mean)
		rmsd_difference = RMSDTools::calcRMS(reference_coords,
				mean_coords,
				this->rmsdData->atomsPerFittingConformation);

		// reference = mean
		RMSDTools::copyArrays(reference_coords,
				mean_coords,
				this->rmsdData->fittingConformationLength);

		return rmsd_difference;
}

void RMSDCalculator::superposition_with_external_reference(double* reference_coords){
	if (! this->rmsdData->hasCalculationCoordinatesSet()){
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
	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			1,
			reference);
	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			-1,
			NULL,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingConformationLength,
			this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->fittingCoordinates);
}

void RMSDCalculator::superposition_with_external_reference_rotating_calc_coords(double* reference){

	// Center fitting coordinates
	double* fitCenters = new double[this->rmsdData->numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
									1,
									reference);

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates,
			fitCenters);

	// Apply a relative translation to center also the coordinates we want to rotate
	// them with the fitting coordinates
	RMSDTools::applyTranslationsToAll(this->rmsdData->atomsPerCalculationConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->calculationCoordinates,
			fitCenters,
			-1);
	delete [] fitCenters;

	this->kernelFunctions->oneVsFollowingFitDiffersCalcWithConfRotation(
			reference,
			NULL,
			-1,
			NULL,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingConformationLength,
			this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->fittingCoordinates,
			this->rmsdData->calculationConformationLength,
			this->rmsdData->atomsPerCalculationConformation,
			this->rmsdData->calculationCoordinates);

}
