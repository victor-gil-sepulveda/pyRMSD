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
			oneVsFollowing(conformation_number,&(rmsd_tmp[offset]));
		}

		for (int i = 0; i < num_of_rmsds; ++i){
			rmsd.push_back(rmsd_tmp[i]);
		}

		delete [] rmsd_tmp;
}

void RMSDomp::oneVsFollowing(int conformation, double *rmsd){
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

void void RMSDomp::iterativeSuperposition(double rmsd_diff_to_stop = 0.001){
	// In the first step, reference is the first conformation
	double* mean_coords = new double[coordinatesPerConformation];

	double* reference = allCoordinates;
	copy(reference, reference+coordinatesPerConformation, mean_coords);

	double last_rmsd = 10000.0; // Big number
	double rmsd_diff = 0.0;
	double	rmsds[numberOfConformations];
	double old_rmsds[numberOfConformations];
	do{
		// Superpose all conformations with the reference conformation
		for (int i = 0; i < numberOfConformations; ++i){
			double* working_conformation = &(allCoordinates[i*coordinatesPerConformation]);
			RMSDTools::superpose(atomsPerConformation, mean_coords, working_conformation);
			rmsds[i] = RMSDTools::calcRMS(mean_coords, working_conformation, atomsPerConformation);
		}
		// Zero mean coordinates
		for (int i  = 0; i <  coordinatesPerConformation; ++i){
			mean_coords[i] = 0;
		}
		// Calculate new mean coords
		for (int i  = 0; i <  numberOfConformations; ++i){
			int conformation_offset = i*coordinatesPerConformation;
			for (int j = 0; j < atomsPerConformation; ++j){
				int atom_offset = 3*j;
				int offset = conformation_offset + atom_offset;
				mean_coords[atom_offset] += allCoordinates[ offset ];
				mean_coords[atom_offset+1] += allCoordinates[ offset + 1];
				mean_coords[atom_offset+2] += allCoordinates[ offset + 2];
			}
		}

	}while(rmsd_diff > rmsd_diff_to_stop);
}

