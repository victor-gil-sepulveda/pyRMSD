#include "ThRMSDSerial.h"
#include "kernel_functions_serial.h"
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
ThRMSDSerial::ThRMSDSerial(int numberOfConformations, int atomsPerConformation, double* coords):
				RMSD(numberOfConformations, atomsPerConformation, coords){
    
	this->tmpRMSDs = new double[numberOfConformations];
	this->kernelFunctions = new ThRMSDSerialKernel;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class destructor. Frees memory.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
ThRMSDSerial::~ThRMSDSerial(){
	delete [] tmpRMSDs;
	delete kernelFunctions;
}

void ThRMSDSerial::_one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double* rmsd){
	// center coords
	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllAtOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	// TODO: evitar copiar el vector de rmsd cambiando la kernel function
	this->kernelFunctions->calcRMSDOfOneVsFollowing(this->allCoordinates, reference_conformation_number, reference_conformation_number + 1,
								numberOfConformations, atomsPerConformation, this->tmpRMSDs);

	// Copy it to the output vector (needs to have the correct size)
	int j = 0;
	for (int i = reference_conformation_number + 1; i < numberOfConformations;++i,++j){
		rmsd[j] = (double) this->tmpRMSDs[i];
	}

	// Move conformations again to their places to avoid coordinate modification
	RMSDTools::applyTranslationsToAll(this->atomsPerConformation, this->numberOfConformations, this->allCoordinates, centers);
	delete [] centers;
}
