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

QTRFITSerialCalculator::~QTRFITSerialCalculator(){
}


void QTRFITSerialCalculator::_one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	// Indeed only numberOfConformations-reference_conformation_number+1 conformations have to be centered...
	// TODO: change this to gain performance
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	dynamic_cast<QTRFITSerialKernel*>(this->getKernelFunctions())->oneVsFollowingFitEqualCalcWithoutConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates);

	// Move then again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, centers);
	delete [] centers;
}


void QTRFITSerialCalculator::_one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd){
	// Center all
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
	delete [] centers;

	dynamic_cast<QTRFITSerialKernel*>(this->getKernelFunctions())->oneVsFollowingFitEqualCalcWithConfRotation(
			reference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates);
	// Coordinates are left centered
}



void QTRFITSerialCalculator::_one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	dynamic_cast<QTRFITSerialKernel*>(this->getKernelFunctions())->oneVsFollowingFitDiffersCalcWithoutConfRotation(
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



void QTRFITSerialCalculator::_one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd){

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllAtOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);
	delete [] fitCenters;
	delete [] calcCenters;

	dynamic_cast<QTRFITSerialKernel*>(this->getKernelFunctions())->oneVsAllFitDiffersCalcWithConfRotation(
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

KernelFunctions* QTRFITSerialCalculator::getKernelFunctions(){
	if(this->kernelFunctions == NULL){
		this->kernelFunctions =  new QTRFITSerialKernel;
	}
	return this->kernelFunctions;
}
