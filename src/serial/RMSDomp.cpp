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

#include <sstream>
#include <fstream>

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

		if(this->allRMSDCoordinates == NULL){
			this->_one_vs_following_fit_equals_calc_coords(conformation, rmsd);
		}
		else{
			this->_one_vs_following_fit_differs_calc_coords(conformation, rmsd);
		}
	}
}

void RMSDomp::_one_vs_following_fit_equals_calc_coords(int conformation, double *rmsd){
	double* reference = &(allCoordinates[conformation*coordinatesPerConformation]);

	double* centers = new double[numberOfConformations*3];
	RMSDTools::centerAllToOrigin(atomsPerConformation, numberOfConformations, allCoordinates, centers);

	#pragma omp parallel for shared(rmsd)
	for (int i = conformation+1; i < numberOfConformations; ++i){

		double* conformation_coords = &(allCoordinates[i*coordinatesPerConformation]);
		double* conformation_coords_tmp = new double[coordinatesPerConformation];
		copy(conformation_coords,conformation_coords+coordinatesPerConformation,conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, conformation_coords_tmp, reference);

		rmsd[i-(conformation+1)] = RMSDTools::calcRMS(reference, conformation_coords_tmp, atomsPerConformation);

		delete [] conformation_coords_tmp;
	}
	delete [] centers;
}

void RMSDomp::_one_vs_following_fit_differs_calc_coords(int conformation, double *rmsd){
	double* fitReference = &(allCoordinates[conformation*coordinatesPerConformation]);
	double* calcReference = &(allRMSDCoordinates[conformation*coordinatesPerRMSDConformation]);

	double* fitCenters = new double[numberOfConformations*3];
	double* calcCenters = new double[numberOfConformations*3];

	RMSDTools::centerAllToOrigin(atomsPerConformation, numberOfConformations, allCoordinates, fitCenters);
	RMSDTools::centerAllToOrigin(atomsPerRMSDConformation, numberOfConformations, allRMSDCoordinates, calcCenters);

	#pragma omp parallel for shared(rmsd)
	for (int i = conformation+1; i < numberOfConformations; ++i){
		int fitCoordsOffset = i*coordinatesPerConformation;
		double* fit_conformation_coords = &(allCoordinates[fitCoordsOffset]);
		double* fit_conformation_coords_tmp = new double[coordinatesPerConformation];
		copy(fit_conformation_coords,fit_conformation_coords+coordinatesPerConformation,fit_conformation_coords_tmp);

		int calcCoordsOffset = i*coordinatesPerRMSDConformation;
		double* calc_conformation_coords = &(allRMSDCoordinates[calcCoordsOffset]);
		double* calc_conformation_coords_tmp = new double[coordinatesPerRMSDConformation];
		copy(calc_conformation_coords, calc_conformation_coords+coordinatesPerRMSDConformation, calc_conformation_coords_tmp);

		RMSDTools::superpose(atomsPerConformation, 	   fit_conformation_coords_tmp,  fitReference,
							 atomsPerRMSDConformation, calc_conformation_coords_tmp, calcReference);

		rmsd[i-(conformation+1)] = RMSDTools::calcRMS(calcReference, calc_conformation_coords_tmp, atomsPerRMSDConformation);

		delete [] fit_conformation_coords_tmp;
		delete [] calc_conformation_coords_tmp;
	}

	delete [] fitCenters;
	delete [] calcCenters;
}



// reference coords is already a copy, in a different memory space than allCoordinates
void RMSDomp::superpositionChangingCoordinates(double* reference, double* rmsds){
	// Superpose all conformations with the reference conformation
	#pragma omp parallel for shared(rmsds)
	for (int i = 0; i < numberOfConformations; ++i){
		double* working_conformation = &(allCoordinates[i*coordinatesPerConformation]);

		double* reference_tmp = new double[coordinatesPerConformation];
		copy(reference, reference+coordinatesPerConformation, reference_tmp);

		RMSDTools::superpose(atomsPerConformation, reference_tmp, working_conformation);

		if (rmsds != NULL){
			rmsds[i] = RMSDTools::calcRMS(reference_tmp, working_conformation, atomsPerConformation);
		}

		delete [] reference_tmp;
	}
}

void RMSDomp::iterativeSuperposition(double rmsd_diff_to_stop = 1e-4){
	// In the first step, reference is the first conformation
	double MAX_ITERATIONS = 200;
	double* reference_coords = new double[coordinatesPerConformation];
	double* mean_coords = new double[coordinatesPerConformation];

	double rmsd_difference = 0.0;
	int current_iteration = 0;

	// reference = coordinates[0]
	RMSDTools::copyArrays(reference_coords, allCoordinates, coordinatesPerConformation);
	do{
		// Superpose all conformations with the reference conformation
		superpositionChangingCoordinates(reference_coords, NULL);

//		for(int i = 0; i < numberOfConformations; ++i ){
//			double* coords = &(allCoordinates[i*coordinatesPerConformation]);
//			save_vector(coords, "coordinates_pyrmsd", coordinatesPerConformation,"Conformation ", i, current_iteration);
//		}

		// Calculate new mean coords, which will be the next reference
		RMSDTools::calculateMeanCoordinates(mean_coords, allCoordinates,
											numberOfConformations, atomsPerConformation);

		// rmsd(reference,mean)
		rmsd_difference = RMSDTools::calcRMS(reference_coords, mean_coords, atomsPerConformation);

		// reference = mean
		RMSDTools::copyArrays(reference_coords, mean_coords, coordinatesPerConformation);

		//cout<< "rmsd diff: "<<current_iteration<<" "<<rmsd_difference<<endl;
		current_iteration++;
	}while(rmsd_difference > rmsd_diff_to_stop and current_iteration < MAX_ITERATIONS);
}
