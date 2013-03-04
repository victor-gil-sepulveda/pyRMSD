#include "ThRMSDSerial.h"
#include "kernel/kernel_functions_serial.h"
#include "../RMSDTools.h"
#include <assert.h>
#include <iostream>

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
///	Class constructor. Allocates memory and copies coordinates with the required precision type.
///
/// \param 	numberOfConformations [In] Total number of conformations in the coordinates array.
///
/// \param 	atomsPerConformation [In] Number of atoms of every conformation.
///
/// \param 	coords [In] Coordinates array with numberOfConformations*atomsPerConformation*3 elements.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
QCPSerialCalculator::QCPSerialCalculator(int numberOfConformations, int atomsPerConformation, double* coords):
				RMSDCalculator(numberOfConformations, atomsPerConformation, coords){
    
	this->kernelFunctions = new ThRMSDSerialKernel;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class destructor. Frees memory.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
QCPSerialCalculator::~QCPSerialCalculator(){
	delete kernelFunctions;
}

void QCPSerialCalculator::_one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double* rmsd){
	// center coords
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	// This kernel function does not modify coordinates.
	dynamic_cast<QCPKernelFunctions*>(this->kernelFunctions)->calcRMSDOfOneVsFollowing(this->allCoordinates,
													reference,
													reference_conformation_number,
													numberOfConformations,
													atomsPerConformation,
													rmsd);

	// Move conformations again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, centers);
	delete [] centers;
}

void QCPSerialCalculator::_one_vs_following_fit_differs_calc_coords(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd){
	// Center coordinates
	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	dynamic_cast<QCPKernelFunctions*>(this->kernelFunctions)->calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoords(
			allCoordinates,
			allRMSDCoordinates,
			fitReference,
			calcReference,
			reference_conformation_number,
			this->numberOfConformations,
			this->atomsPerConformation,
			this->atomsPerRMSDConformation,
			rmsd);

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, fitCenters);
	RMSDTools::applyTranslationsToAll(this->atomsPerRMSDConformation, this->numberOfConformations, this->allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;
}


void QCPSerialCalculator::_one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
	delete [] centers;
	dynamic_cast<QCPKernelFunctions*>(this->kernelFunctions)->calcRMSDOfOneVsFollowingModifyingCoordinates(this->allCoordinates,
																		reference,
																		reference_conformation_number,
																		numberOfConformations,
																		atomsPerConformation,
																		rmsd);
}

void QCPSerialCalculator::_one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference,
		double* calcReference, int reference_conformation_number, double *rmsd){
		double* fitCenters = new double[numberOfConformations*3];
		double* calcCenters = new double[numberOfConformations*3];
		RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
		RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);
		delete [] fitCenters;
		delete [] calcCenters;

		dynamic_cast<QCPKernelFunctions*>(this->kernelFunctions)->calcRMSDOfOneVsFollowingWithDifferentFitAndCalcCoordsModifyingCoordinates(
				allCoordinates,
				allRMSDCoordinates,
				fitReference,
				calcReference,
				reference_conformation_number,
				this->numberOfConformations,
				this->atomsPerConformation,
				this->atomsPerRMSDConformation,
				rmsd);

}
