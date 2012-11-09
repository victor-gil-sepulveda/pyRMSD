#include "ThRMSDSerial.h"
#include "kernel_functions_serial.h"
#include <assert.h>
#include <iostream>

using namespace std;

#define floating_point_type double


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

ThRMSDSerial::~ThRMSDSerial(){
    delete [] this->convertedCoords;
    delete [] this->tmpRMSDs;
}

void ThRMSDSerial::oneVsTheOthers(int conformation, double* rmsd_result) {
	if (conformation < numberOfConformations){
		//cout<<conformation<<" vs the others"<<endl;
		
		ThRMSDSerialKernel::calcRMSDOfOneVsOthers(convertedCoords, conformation, conformation + 1,
	    			numberOfConformations, atomsPerConformation, this->tmpRMSDs);
	    
	    // Do the copy to the output vector (need to have the correct size)
		int j = 0;
		for (int i = conformation + 1; i < numberOfConformations;++i,++j){
			rmsd_result[j] = (double) this->tmpRMSDs[i];
		}
	}
}

void ThRMSDSerial::calculateRMSDCondensedMatrix(vector<double>& rmsd){
	for(int conformation_number = 0;conformation_number < numberOfConformations; ++conformation_number){ 
		int number_of_rmsds = numberOfConformations - conformation_number-1;
		double* rmsd_tmp = new double[number_of_rmsds];
		oneVsTheOthers(conformation_number,rmsd_tmp);
		for (int i = 0; i < number_of_rmsds; ++i){
			rmsd.push_back(rmsd_tmp[i]);
		}
		delete [] rmsd_tmp;
	}
}

double* ThRMSDSerial::getCoordinates(){
	return convertedCoords;
}
