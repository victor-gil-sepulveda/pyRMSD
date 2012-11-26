#include "ThRMSDSerial.h"
#include "kernel_functions_serial.h"
#include <assert.h>
#include <iostream>

using namespace std;

#define floating_point_type double

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
    
	// Allocate space for temporary rmsds
	this->tmpRMSDs = new floating_point_type[numberOfConformations];
	
	int total_num_of_coords = this->numberOfConformations*this->atomsPerConformation*3;

	// Allocate space and convert coords
	this->convertedCoords = new floating_point_type[total_num_of_coords];
	for (int i =0; i<total_num_of_coords;++i){
		this->convertedCoords[i] = (floating_point_type) coords[i];
	}

	ThRMSDSerialKernel::centerCoordsOfAllConformations(numberOfConformations,atomsPerConformation, convertedCoords);
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class destructor. Frees memory.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
ThRMSDSerial::~ThRMSDSerial(){
    delete [] this->convertedCoords;
    delete [] this->tmpRMSDs;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Implements a 'i' vs [i+1,M] rmsd conformation calculation, where M is the number of conformations.
///
/// \param 	conformation [In] Total number of conformations in the coordinates array.
///
/// \param 	rmsd_result [In/Out] Array where the rmsds will be stored (it needs to have the correct size!).
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDSerial::oneVsFollowing(int conformation, double* rmsd_result) {
	if (conformation < numberOfConformations){
		//cout<<conformation<<" vs the others"<<endl;
		
		ThRMSDSerialKernel::calcRMSDOfOneVsOthers(convertedCoords, conformation, conformation + 1,
	    			numberOfConformations, atomsPerConformation, this->tmpRMSDs);
	    
	    // Do the copy to the output vector (needs to have the correct size)
		int j = 0;
		for (int i = conformation + 1; i < numberOfConformations;++i,++j){
			rmsd_result[j] = (double) this->tmpRMSDs[i];
		}
	}
}


///////////////////////////////////////////////////////////////
/// \remarks
///	Implements the calculation of the upper triangular matrix in row major format of the rmsd distance between all
/// conformations stored in the coordinates.
///
/// \param 	rmsd [In/Out] Vector where the resulting RMSDs will be stored.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDSerial::calculateRMSDCondensedMatrix(vector<double>& rmsd){
	for(int conformation_number = 0;conformation_number < numberOfConformations; ++conformation_number){ 
		int number_of_rmsds = numberOfConformations - conformation_number-1;
		double* rmsd_tmp = new double[number_of_rmsds];
		oneVsFollowing(conformation_number,rmsd_tmp);
		for (int i = 0; i < number_of_rmsds; ++i){
			rmsd.push_back(rmsd_tmp[i]);
		}
		delete [] rmsd_tmp;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Returns the class representation of coordinates (for testing purposes, only works with double version).
///
/// \return	The pointer to the array containing converted coordinates.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
double* ThRMSDSerial::getCoordinates(){
	return convertedCoords;
}
