/*
 * RMSDomp.cpp
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#include "RMSDomp.h"
#include "RMSDTools.h"
#include <omp.h>
#include <vector>
#include <iostream>
#include <algorithm>
using namespace std;

RMSDomp::RMSDomp(int numberOfConformations, int atomsPerConformation, double *allCoordinates):
		RMSD(numberOfConformations, atomsPerConformation,allCoordinates){
}

RMSDomp::~RMSDomp() {
}

void RMSDomp::calculateRMSDCondensedMatrix(std::vector<double> & rmsd){
	int num_of_rmsds = numberOfConformations * (numberOfConformations-1.) /2.;
		double* rmsd_tmp =  new double[num_of_rmsds];

		for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){
			int offset = conformation_number*(numberOfConformations-1)- (((conformation_number-1)*conformation_number)/2) ;
			oneVsTheOthers(conformation_number,&(rmsd_tmp[offset]));
		}

		for (int i = 0; i < num_of_rmsds; ++i){
			rmsd.push_back(rmsd_tmp[i]);
		}

		delete [] rmsd_tmp;
}

void RMSDomp::oneVsTheOthers(int conformation, double *rmsd){
	if(conformation>=numberOfConformations){
		cout<<"Error, this conformation doesn't exist ("<<conformation<<")"<<endl;
	}
	else{
		int coordsOffset = conformation*coordinatesPerConformation;
		double* reference = &(allCoordinates[coordsOffset]);


		#pragma omp parallel for shared(rmsd)
		for (int i = conformation+1; i < numberOfConformations; ++i){
			int rmsd_index = i-(conformation+1);

			coordsOffset = i*coordinatesPerConformation;

			double* conformation_coords = &(allCoordinates[coordsOffset]);

			double* reference_tmp = new double[coordinatesPerConformation];
			double* conformation_coords_tmp = new double[coordinatesPerConformation];

			copy(reference,reference+coordinatesPerConformation,reference_tmp);
			copy(conformation_coords,conformation_coords+coordinatesPerConformation,conformation_coords_tmp);

			RMSDTools::superpose(atomsPerConformation,reference_tmp,conformation_coords_tmp);
			double rmsd_val = RMSDTools::calcRMS(reference_tmp,conformation_coords_tmp,atomsPerConformation);

			delete [] reference_tmp;
			delete [] conformation_coords_tmp;

			rmsd[rmsd_index] = rmsd_val;
		}
	}
}


