#include "ThRMSDSerialOmp.h"
#include "ThRMSDSerial.h"
#include "kernel_functions_omp.h"
#include "kernel_functions_serial.h"
#include <assert.h>
#include <iostream>
#include <omp.h>

using namespace std;

#define floating_point_type double


ThRMSDSerialOmp::ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation,
		double* coords, int omp_threads):
				ThRMSDSerial(numberOfConformations, atomsPerConformation, coords){
    this->omp_threads = omp_threads;
	omp_set_num_threads(this->omp_threads);
	cout<<"Setting omp threads to "<<this->omp_threads<<endl;
}

void ThRMSDSerialOmp::oneVsFollowing(int conformation, double* rmsd_result) {
	if (conformation < numberOfConformations){
		
		ThRMSDSerialOmpKernel::calcRMSDOfOneVsFollowing(convertedCoords, conformation, conformation + 1,
	    			numberOfConformations, atomsPerConformation, this->tmpRMSDs, this->omp_threads);
	    
	    // Do the copy to the output vector (needs to have the correct size)
		int j = 0;
		for (int i = conformation + 1; i < numberOfConformations;++i,++j){
			rmsd_result[j] = (double) this->tmpRMSDs[i];
		}
	}
}
