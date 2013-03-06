#include <cmath>
#include <iostream>
#include "QCPCUDAKernel.h"
#include "../RMSDTools.h"
#include "kernel_functions_cuda.h"

using namespace std;

///////////////////////////////////////////////////////////////
/// \remarks
///	Convenience function for CUDA error handling. It captures the error, writes a user-based message,
/// and exits the program.
///
/// \param 	message [In] Message to print if something went wrong (usually calling CUDA function name).
///
/// \param 	error_code [In] CUDA error code (It may be 0 if everything is OK).
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
inline void checkCudaError(char* message, cudaError_t error_code){
	if (error_code != 0){
		std::cout<<"Error in "<<message<<" . Error code: "<<error_code<<". Exiting..."<<std::flush<<std::endl;
		exit(-1);
	}
}

void QCPCUDAKernel::init(
		double* coordinates, 
		int atomsPerConformation,
		int coordinatesPerConformation, 
		int numberOfConformations){

	cudaDeviceProp props;
	int device;
	checkCudaError("cudaGetDevice",cudaGetDevice(&device));
	checkCudaError("cudaGetDeviceProperties",cudaGetDeviceProperties(&props, device));

	if(!props.deviceOverlap){
		cout<<"There are no multiple streams, so this implementation may not work properly :S"<<endl;
	}

	// Convert to our floating point type
	int totalNumberOfCoordinates = coordinatesPerConformation*numberOfConformations;
	this->tmpHostCoords = new floating_point_type[totalNumberOfCoordinates];
	for(int i = 0; i < totalNumberOfCoordinates; ++i){
		this->tmpHostCoords[i] = (floating_point_type) coordinates[i];
	}

	// Allocate space for temporary rmsds
	this->tmpHostRMSDs = new floating_point_type[numberOfConformations];

	// GPU Data allocation for input
	checkCudaError(" Malloc Device Coords. ", cudaMalloc((void **) &deviceCoords, totalNumberOfCoordinates * sizeof(floating_point_type)));
	checkCudaError(" Copying Coords to Device. ", cudaMemcpy(deviceCoords, tmpHostCoords, totalNumberOfCoordinates * sizeof(floating_point_type),
			cudaMemcpyHostToDevice));

	// GPU Data allocation for output
	checkCudaError(" Malloc RMSDs ", cudaMalloc((void **) &deviceRMSDs, numberOfConformations * sizeof(floating_point_type)));

	centerCoordsOfAllConformations<<<this->blocks_per_grid, this->threads_per_block>>>(
			numberOfConformations,
			atomsPerConformation,
			deviceCoords);
}

QCPCUDAKernel::~QCPCUDAKernel(){
	checkCudaError(" deviceCoords cudaFree ", cudaFree(this->deviceCoords));
	checkCudaError(" deviceRMSDs cudaFree ", cudaFree(this->deviceRMSDs));

	if (tmpHostCoords != NULL){
		delete [] tmpHostCoords;
	}

	if (tmpHostRMSDs != NULL){
		delete [] tmpHostRMSDs;
	}
}

void QCPCUDAKernel::changeCalculationCoords(double*){

}

void QCPCUDAKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates){

	calcRMSDOfOneVsFollowing<<<this->blocks_per_grid, this->threads_per_block>>>(
			 this->deviceCoords,
			 reference_conformation_number,
			 reference_conformation_number+1,
			 numberOfConformations,
			 atomsPerConformation,
			 coordinatesPerConformation,
			 this->deviceRMSDs);

	checkCudaError(" Getting RMSDs from Device. ", cudaMemcpy(this->tmpHostRMSDs, this->deviceRMSDs,
												   numberOfConformations * sizeof(floating_point_type),
												   cudaMemcpyDeviceToHost));

	// Do the copy to the output vector (needs to have the correct size)
	int j = 0;
	for (int i = reference_conformation_number + 1; i < numberOfConformations;++i,++j){
		rmsd[j] = (double) this->tmpHostRMSDs[i];
	}
}

void QCPCUDAKernel::oneVsFollowingFitDiffersCalcWithoutConfRotation(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates,
		int coordinatesPerRMSDConformation,
		int atomsPerRMSDConformation,
		double *allRMSDCoordinates){


}


void QCPCUDAKernel::oneVsFollowingFitEqualCalcWithConfRotation(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates){

}

void QCPCUDAKernel::oneVsAllFitDiffersCalcWithConfRotation(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates,
		int coordinatesPerRMSDConformation,
		int atomsPerRMSDConformation,
		double *allRMSDCoordinates){

}

