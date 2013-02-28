#include "RMSDSerial.h"
#include "RMSDTools.h"
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

RMSDSerial::RMSDSerial(int numberOfConformations, int atomsPerConformation, double* allCoordinates):
		RMSD(numberOfConformations, atomsPerConformation,allCoordinates){
}

RMSDSerial::~RMSDSerial(){}
 
void RMSDSerial::oneVsFollowing(int conformation, double* rmsd){
	if(conformation>=numberOfConformations){
		cout<<"Error, this conformation doesn't exist ("<<conformation<<")"<<endl;
	}
	else{
		int coordsOffset = conformation*coordinatesPerConformation;
		double* reference = &(allCoordinates[coordsOffset]);
		
		double* centers = new double[numberOfConformations*3];
		RMSDTools::centerAllToOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);
		delete [] centers;

		int rmsd_index = 0;
		for (int i = conformation+1; i < numberOfConformations; ++i){
			coordsOffset += coordinatesPerConformation;
			double* conformation_coords = &(allCoordinates[coordsOffset]);
			RMSDTools::superpose(atomsPerConformation,reference,conformation_coords);
			double rmsd_val = RMSDTools::calcRMS(reference,conformation_coords,atomsPerConformation);
			rmsd[rmsd_index] = rmsd_val;
			rmsd_index++;
		}

	}
}

void RMSDSerial::calculateRMSDCondensedMatrix(vector<double>& rmsd){ 
	for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){ 
		//cout<<"Calculating "<<conformation_number<<"th conformation ("<<conformation_number*100./numberOfConformations<<"%)"<<endl<<flush;
		int number_of_rmsds = numberOfConformations-conformation_number-1;
		double* rmsd_tmp = new double[number_of_rmsds];
    	oneVsFollowing(conformation_number,rmsd_tmp);
		for (int i = 0; i < number_of_rmsds; ++i){
			rmsd.push_back(rmsd_tmp[i]);
		}
		delete [] rmsd_tmp;
	}
}
