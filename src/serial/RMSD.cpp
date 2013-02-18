#include "RMSD.h"
#include <iostream>
using namespace std;

RMSD::RMSD(int numberOfConformations, int atomsPerConformation, double* allCoordinates){
	this->numberOfConformations = numberOfConformations;
	this->atomsPerConformation = atomsPerConformation;
	this->allCoordinates = allCoordinates;
	this->coordinatesPerConformation = atomsPerConformation*3;
	this->allRMSDCoordinates = NULL;
}

/*
 * Sets a different set of coordinates for RMSD calculation and fit.
 */
void RMSD::setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates){
	this->atomsPerRMSDConformation = atomsPerRMSDConformation;
	this->coordinatesPerRMSDConformation = atomsPerRMSDConformation*3;
	this->allRMSDCoordinates = allRMSDCoordinates;
}

void RMSD::iterativeSuperposition(double rmsd_diff_to_stop = 1e-4){
	cout<<"Iterative superposition not implemented yet for this algorithm."<<flush<<endl;
}
