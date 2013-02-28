#include "RMSD.h"
#include <iostream>
#include <cstdlib>
using namespace std;

RMSD::RMSD(int numberOfConformations, int atomsPerConformation, double* allCoordinates){
	this->numberOfConformations = numberOfConformations;
	this->atomsPerConformation = atomsPerConformation;
	this->allCoordinates = allCoordinates;
	this->coordinatesPerConformation = atomsPerConformation*3;
	this->allRMSDCoordinates = NULL;
	this->coordinatesPerRMSDConformation = 0;
	this->atomsPerRMSDConformation = 0;
	this->modifyFittingCoordinates = false;
}

/*
 * Sets a different set of coordinates for RMSD calculation and fit.
 */
void RMSD::setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates){
	this->atomsPerRMSDConformation = atomsPerRMSDConformation;
	this->coordinatesPerRMSDConformation = atomsPerRMSDConformation*3;
	this->allRMSDCoordinates = allRMSDCoordinates;
}

/*
 *
 * Iterative superoposition algorithm
 *
 */
void RMSD::iterativeSuperposition(double rmsd_diff_to_stop = 1e-4){
	cout<<"[Error RMSD::Core::iterativeSuperposition]Iterative superposition not implemented for this calculator."<<flush<<endl;
	exit(-1);
}
