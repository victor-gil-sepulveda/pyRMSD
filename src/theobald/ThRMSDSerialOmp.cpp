#include "ThRMSDSerialOmp.h"
#include "ThRMSDSerial.h"
#include "kernel_functions_omp.h"
#include "kernel_functions_serial.h"
#include <assert.h>
#include <iostream>

using namespace std;

#define floating_point_type double


ThRMSDSerialOmp::ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation, double* coords):
				ThRMSDSerial(numberOfConformations, atomsPerConformation, coords){
    
}

void ThRMSDSerialOmp::oneVsFollowing(int conformation, double* rmsd_result) {
	if (conformation < numberOfConformations){
		//cout<<conformation<<" vs the others"<<endl;
		
		ThRMSDSerialOmpKernel::calcRMSDOfOneVsOthers(convertedCoords, conformation, conformation + 1,
	    			numberOfConformations, atomsPerConformation, this->tmpRMSDs);
	    
	    // Do the copy to the output vector (need to have the correct size)
		int j = 0;
		for (int i = conformation + 1; i < numberOfConformations;++i,++j){
			rmsd_result[j] = (double) this->tmpRMSDs[i];
		}
	}
}
