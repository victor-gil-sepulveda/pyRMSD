#include "RMSDCalculator.h"
#include "RMSDTools.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
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
		this->kernelFunctions->setCalculationCoords(rmsdData);
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
			this->_one_vs_following_fit_equals_calc_coords(
					reference_conformation,
					reference_conformation_index,
					rmsd);

		}
		else{
			double* fit_reference_conformation = this->rmsdData->getFittingConformationAt(reference_conformation_index);
			double* calc_reference_conformation = this->rmsdData->getCalculationConformationAt(reference_conformation_index);

			this->_one_vs_following_fit_differs_calc_coords(
					fit_reference_conformation,
					calc_reference_conformation,
					reference_conformation_index,
					rmsd);
			if(this->rmsdData->hasSymmetryGroups()){
				this->kernelFunctions->handleSymmetriesWithCalcCoords(
												reference_conformation_index,
												rmsd,
												this->rmsdData);
			}

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
		calculate_rmsd_condensed_matrix_with_fitting_coordinates(rmsd);
	}
	else{
		calculate_rmsd_condensed_matrix_with_fitting_and_calculation_coordinates(rmsd);
	}

	// Post processing (for CUDA based, for instance)
	this->kernelFunctions->matrixEnd(num_of_rmsds,rmsd);
}

/**
 * Does the actual calculation of the matrix, doing the RMSD calculation using the
 * same coordinates that were used for fitting.
 *
 * \param rmsd Linear vector that will hold the rmsd values of the matrix. The Nth-1 elements will
 * form the first row, the next N-2 the second row and so on.
 */
void RMSDCalculator::calculate_rmsd_condensed_matrix_with_fitting_coordinates(vector<double>& rmsd){
	
	this->kernelFunctions->centerAllAtCOM(this->rmsdData);
	
	this->kernelFunctions->matrixInit(rmsdData);

	for(int reference_index = 0; reference_index < this->rmsdData->numberOfConformations - 1; ++reference_index){
			int offset = reference_index*(this->rmsdData->numberOfConformations-1)- (((reference_index-1)*reference_index)/2) ;
			double* reference_conformation = this->rmsdData->getFittingConformationAt(reference_index);
			this->kernelFunctions->matrixOneVsFollowingFitEqualCalc(
					reference_conformation,
					reference_index,
					&(rmsd[offset]),
					rmsdData);
	}
}

/**
 * Exactly the same than 'calculateRMSDCondensedMatrixWithFittingCoordinates' but this time uses
 * one coordinate set for fitting and other to calculate the RMSD.
 *
 * \param rmsd Linear vector that will hold the rmsd values of the matrix. The Nth-1 elements will
 * form the first row, the next N-2 the second row and so on.
 */
void RMSDCalculator::calculate_rmsd_condensed_matrix_with_fitting_and_calculation_coordinates(vector<double>& rmsd){

	this->kernelFunctions->centerAllAtFittingCOM(this->rmsdData);
	
	this->kernelFunctions->matrixInit(rmsdData);

	for(int reference_index = 0; reference_index<this->rmsdData->numberOfConformations; ++reference_index){
		int offset = reference_index*(this->rmsdData->numberOfConformations-1)- (((reference_index-1)*reference_index)/2) ;
		double* fit_reference_conformation = this->rmsdData->getFittingConformationAt(reference_index);
		double* calc_reference_conformation = this->rmsdData->getCalculationConformationAt(reference_index);
		this->kernelFunctions->matrixOneVsFollowingFitDiffersCalc(
															fit_reference_conformation,
															calc_reference_conformation,
															reference_index,
															&(rmsd[offset]),
															rmsdData);

		if(this->rmsdData->hasSymmetryGroups()){
			this->kernelFunctions->handleSymmetriesWithCalcCoords(
											reference_index,
											&(rmsd[offset]),
											this->rmsdData);
		}
	}

}

/**
 *	Prepares coordinates and launches the kernel function that will do the i vs i+1:N superpostion +
 *	rmsd calculation operation where the RMSD calculation will be performed with the fitting coordinates.
 *	Modifies all input coordinates.
 *
 * \param reference A pointer to the reference conformation.
 * \param reference_index The place that the reference conformation has in the conformations sequence
 * (the fitting coordinates array).
 * \param rmsd An array with enough room to store the calculated rmsd values.
 */
void RMSDCalculator::_one_vs_following_fit_equals_calc_coords(
		double* reference,
		int reference_index,
		double *rmsd){

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcCoords(
			reference,
			reference_index,
			rmsd,
			rmsdData);
}

/**
 *	Prepares coordinates and launches the kernel function that will do the i vs i+1:N superpostion +
 *	rmsd calculation operation. The RMSD calculation will be performed over a set which is different
 *	from the fitting coordinates set. Modifies all input coordinates.
 *
 * \param reference A pointer to the reference conformation.
 * \param reference_index The place that the reference conformation has in the conformations sequence
 * (the fitting coordinates array).
 * \param rmsd An array with enough room to store the calculated rmsd values.
 */
void RMSDCalculator::_one_vs_following_fit_differs_calc_coords(
		double* fitReference,
		double* calcReference,
		int reference_index,
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

	this->kernelFunctions->oneVsFollowingFitDiffersCalcCoords(
			fitReference,
			calcReference,
			reference_index,
			rmsd,
			rmsdData);
}

/**
 *	Performs an iterative superposition algorithm (it tries to be equivalent to the one implemented in
 *	prody @ http://www.csb.pitt.edu/prody/). If the calculation coordinates set has been defined, it will
 *	also rotate them so that one can define a set of coordinates for fitting, and a set of coordinates
 *	for querying later. Input coordinates are modified.
 *
 *	Prody scheme is:
 *
 *	 Ensemble @ .../prody/ensemble/ensemble.py
 *	 PDBEnsemble @ .../prody/ensemble/pdbensemble.py
 *
 *	 Steps:
 *	 1: PDBEnsemble::iterpose -> confs_tmp = confs
 *	 2: PDBEnsemble::iterpose -> Ensemble::iterpose(confs_tmp)
 *	 Iterative
 *	 		3: PDBEnsemble::_superpose()
 *	 4: Ensemble::superpose()
 *	 5: PDBEnsemble::_superpose(trans=True)
 *
 * \param rmsd_diff_to_stop Once the rmsd difference is lower than this value the algorithm has
 * converged.
 * \param iteration_rmsd Correctly sized vector where the rmsd difference value of each step will
 * be stored.
 */
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

/**
 * Performs one step of the iterative superposition method.
 *
 * \param reference_coords A vector with the reference conformation. It is the first conformation in the
 * first step, and mean conformation in next steps.
 * \param mean_coords A buffer vector to calculate the mean conformation.
 */
double RMSDCalculator::iterative_superposition_step(double* reference_coords, double* mean_coords){
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


/**
 * Helper function to choose the superposition method depending on the definition or not of the
 * extra calculation set.
 *
 * \param reference Is the reference conformation to be used in superposition (mean conformation
 * in this case).
 */
void RMSDCalculator::superposition_with_external_reference(double* reference){
	if (! this->rmsdData->hasCalculationCoordinatesSet()){
		superposition_with_external_reference_without_calc_coords(reference);
	}
	else{
		superposition_with_external_reference_rotating_calc_coords(reference);
	}
}

/**
 * Superposes all coordinates to the reference conformation. In this case reference coordinates can be
 * placed into a different memory space than the working coordinates set.
 *
 * \param reference Is the reference conformation to be used in superposition (mean conformation
 * in this case).
 */
void RMSDCalculator::superposition_with_external_reference_without_calc_coords(double* reference){

	// Recentering coordinates each step gives us am easier way of comparing with prody,
	// and removes any translation artifact
	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			1,
			reference);

	RMSDTools::centerAllAtOrigin(this->rmsdData->atomsPerFittingConformation,
			this->rmsdData->numberOfConformations,
			this->rmsdData->fittingCoordinates);

	this->kernelFunctions->oneVsFollowingFitEqualCalcCoords(
			reference,
			-1,
			NULL,
			rmsdData);
}

/**
 * Exactly the same that 'superposition_with_external_reference_without_calc_coords', but in this case
 * it will also modify the calculation coordinates set.
 *
 * \param reference Is the reference conformation to be used in superposition ( mean conformation
 * in this case).
 */
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

	this->kernelFunctions->oneVsFollowingFitDiffersCalcCoords(
			reference,
			NULL,
			-1,
			NULL,
			rmsdData);
}
