#include "RMSDCUDA.cuh"
#include "kernel_rmsd.cuh"
#include <assert.h>
#include <iostream>
using namespace std;

inline void checkCudaError(char* id, cudaError_t error_code){
	if (error_code != 0){
		cout<<"Error in "<<id<<" . Error code: "<<error_code<<". Exiting..."<<flush<<endl;
		exit(-1);
	}
}


RMSDCuda::RMSDCuda(int numberOfConformations, int atomsPerConformation, double* coords, int threads_per_block):
				RMSD(numberOfConformations, atomsPerConformation, coords){
    
    //////////////////////////////////////////////////////////
	//// Preparing coordinates for increased coalescence
	//////////////////////////////////////////////////////////
	int total_num_of_coords = this->numberOfConformations*this->atomsPerConformation*3;
	this->coalescedCoords = new float[total_num_of_coords];
	for (int confId =0; confId < numberOfConformations; ++confId){
	    for (int atomId = 0; atomId < atomsPerConformation; ++atomId){
	    	int offset = confId * atomsPerConformation*3 + atomId*3;
	    	this->coalescedCoords[atomId*numberOfConformations*3+0*numberOfConformations+confId] = (float) coords[offset];
	    	this->coalescedCoords[atomId*numberOfConformations*3+1*numberOfConformations+confId] = (float) coords[offset+1];
	    	this->coalescedCoords[atomId*numberOfConformations*3+2*numberOfConformations+confId] = (float) coords[offset+2];
	    }
	 }
	
	// Allocate space for temporary rmsds
	this->tmpHostRMSDs = new float[numberOfConformations];
    
	//////////////////////////////////////////////////////////
	//// Set up GPU
	//////////////////////////////////////////////////////////
	threadsPerBlock = threads_per_block;
	numBlocks = ceil( (1.*numberOfConformations) / threadsPerBlock );

	// GPU Data allocation for input
	checkCudaError(" Malloc Device Coords. ", cudaMalloc((void **) &deviceCoords, total_num_of_coords*sizeof(float)));
	checkCudaError(" Copying Coords to Device. ", cudaMemcpy(deviceCoords,coalescedCoords, total_num_of_coords*sizeof(float), 
			cudaMemcpyHostToDevice));
    
    // GPU Data allocation for output
    checkCudaError(" Malloc Rmsds ", cudaMalloc((void **) &deviceRMSDs, numberOfConformations * sizeof(float)));
    
    //////////////////////////////////////////////////////////
	//// And do some pre-processing
	//////////////////////////////////////////////////////////
    k_center_conformers<<<numBlocks, threadsPerBlock>>>(numberOfConformations,atomsPerConformation,deviceCoords);
}

RMSDCuda::~RMSDCuda(){
    cudaFree(this->deviceCoords);
    cudaFree(this->deviceRMSDs);
    delete [] this->tmpHostRMSDs;
    delete [] this->coalescedCoords;
}

void RMSDCuda::oneVsTheOthers(int conformation, double* rmsd_result) {
    k_all_against_one<<<numBlocks,threadsPerBlock>>>(atomsPerConformation, numberOfConformations, conformation, this->deviceCoords, this->deviceRMSDs);
    
    checkCudaError(" Getting Coords from Device. ", cudaMemcpy(this->tmpHostRMSDs, this->deviceRMSDs, numberOfConformations * sizeof(float), cudaMemcpyDeviceToHost));
    
    // Do the copy to the output vector (need to have the correct size)
	int j = 0;
	for (int i = conformation + 1; i < numberOfConformations;++i,++j){
		rmsd_result[j] = (double) this->tmpHostRMSDs[i];
	}
}

void RMSDCuda::calculateRMSDCondensedMatrix(vector<double>& rmsd){ 
	for(int conformation_number = 0;conformation_number<numberOfConformations;++conformation_number){ 
		int number_of_rmsds = numberOfConformations - conformation_number-1;
		double* rmsd_tmp = new double[number_of_rmsds];
    	oneVsTheOthers(conformation_number,rmsd_tmp);
		for (int i = 0; i < number_of_rmsds; ++i){
			rmsd.push_back((double)rmsd_tmp[i]);
		}
		delete [] rmsd_tmp;
	}
}
