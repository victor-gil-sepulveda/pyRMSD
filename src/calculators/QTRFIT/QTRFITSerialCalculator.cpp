#include "QTRFITSerialCalculator.h"
#include "../RMSDTools.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "kernel/QTRFITSerialKernel.h"
using namespace std;

QTRFITSerialCalculator::QTRFITSerialCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates):
		RMSDCalculator(numberOfConformations, atomsPerConformation,allCoordinates){
}

QTRFITSerialCalculator::~QTRFITSerialCalculator(){}


void QTRFITSerialCalculator::_one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	// Indeed only numberOfConformations-reference_conformation_number+1 conformations have to be centered...
	// TODO: change this to gain performance
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	oneVsFollowingFitEqualCalcWithoutConfRotation(reference_conformation_number,
			reference, rmsd);

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, centers);
	delete [] centers;
}


void QTRFITSerialCalculator::_one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
	delete [] centers;

	oneVsFollowingFitEqualCalcWithConfRotation(reference_conformation_number,
			reference, rmsd);
	// Coordinates are left centered
}



void QTRFITSerialCalculator::_one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	oneVsFollowingFitDiffersCalcWithoutConfRotation(
			reference_conformation_number, fitReference, rmsd, calcReference);

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, fitCenters);
	RMSDTools::applyTranslationsToAll(this->atomsPerRMSDConformation, this->numberOfConformations, this->allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;
}



void QTRFITSerialCalculator::_one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;

	oneVsAllFitDiffersCalcWithConfRotation(reference_conformation_number,
			fitReference, rmsd, calcReference);
}

KernelFunctions* QTRFITSerialCalculator::getKernelFunctions(){
	return new QTRFITSerialKernel;
}
