#include <cmath>
#include <iostream>
#include "kernel_functions_cuda.h"
#include "QCPCUDAKernel.h"
#include "../RMSDTools.h"

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

	int device;
	cudaDeviceProp props;
	
	checkCudaError("cudaGetDevice",
			cudaGetDevice(
					&device));
	
	checkCudaError("cudaGetDeviceProperties",
			cudaGetDeviceProperties(
					&props, 
					device));

	if(!props.deviceOverlap){
		cout<<"There are no multiple streams, so this implementation may not work properly :S"<<endl;
	}

	int totalNumberOfCoordinates = coordinatesPerConformation*numberOfConformations;

	#ifdef CUDA_PRECISION_SINGLE
		// Allocate space for temporary coords
		this->tmpHostCoords = new float[totalNumberOfCoordinates];
		this->tmpHostRMSDs = new float[numberOfConformations];
	#else
		this->tmpHostCoords = NULL;
		this->tmpHostRMSDs = NULL;//new double[numberOfConformations];
	#endif
		
	// GPU Data allocation for input
	checkCudaError(" Malloc Device Coords. ", 
			cudaMalloc(
					(void **) &deviceCoords, 
					totalNumberOfCoordinates * sizeof(floating_point_type)));
	

	// GPU Data allocation for output
	checkCudaError(" Malloc RMSDs ", 
			cudaMalloc(
					(void **) &deviceRMSDs, 
					numberOfConformations * sizeof(floating_point_type)));

	/*centerCoordsOfAllConformations<<<this->blocks_per_grid, this->threads_per_block>>>(
			numberOfConformations,
			atomsPerConformation,
			deviceCoords);*/
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

void QCPCUDAKernel::updateDeviceCoordinates(
		int numberOfConformations,
		int coordinatesPerConformation,
		double * allCoordinates
		){
	
	int totalNumberOfCoordinates = coordinatesPerConformation*numberOfConformations;
	
	#ifdef CUDA_PRECISION_SINGLE
		// Convert to our floating point type
		
		for(int i = 0; i < totalNumberOfCoordinates; ++i){
			this->tmpHostCoords[i] = (float) allCoordinates[i];
		}
		checkCudaError(" Copying Coords to Device (single). ", 
					cudaMemcpy(
							deviceCoords, 
							tmpHostCoords, 
							totalNumberOfCoordinates * sizeof(float),
							cudaMemcpyHostToDevice));
	#else	
		checkCudaError(" Copying Coords to Device (double). ", 
				cudaMemcpy(
						deviceCoords, 
						allCoordinates, 
						totalNumberOfCoordinates * sizeof(double),
						cudaMemcpyHostToDevice));
	#endif
}

void QCPCUDAKernel::updateHostRMSDs(
				int numberOfConformations,
				int reference_conformation_number,
				double* rmsd){
	
	#ifdef CUDA_PRECISION_SINGLE
		// Get RMSDs from device
		checkCudaError(" Getting RMSDs from Device (single). ", 
				cudaMemcpy(
						this->tmpHostRMSDs,
						this->deviceRMSDs,
						(numberOfConformations - reference_conformation_number-1) * sizeof(float),
						cudaMemcpyDeviceToHost));
		
		// Apply conversion
		for (int i = 0 ; i < numberOfConformations - reference_conformation_number-1;++i){
			rmsd[i] = (double) this->tmpHostRMSDs[i];
		}
		
	#else	
		checkCudaError(" Getting RMSDs from Device (double). ", 
				cudaMemcpy(
						rmsd,
						this->deviceRMSDs,
						(numberOfConformations - reference_conformation_number) * sizeof(double),
						cudaMemcpyDeviceToHost));
	#endif

}

void QCPCUDAKernel::oneVsFollowingFitEqualCalcWithoutConfRotation(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates){
	
	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			numberOfConformations,
			coordinatesPerConformation,
			allCoordinates);

	// Do the calculations
	calcRMSDOfOneVsFollowing<<<this->blocks_per_grid, this->threads_per_block>>>(
			 this->deviceCoords,
			 reference_conformation_number,
			 reference_conformation_number+1,
			 numberOfConformations,
			 atomsPerConformation,
			 coordinatesPerConformation,
			 this->deviceRMSDs);

	// Get RMSDs
	updateHostRMSDs(
				numberOfConformations,
				reference_conformation_number,
				rmsd);
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
		/*
		// Update the reference (it can be from another coordinate set)
		updateReferenceCoordinates(
				numberOfConformations,
				coordinatesPerConformation,
				allCoordinates);
	
		// Put the centered coordinates on the device
		updateDeviceCoordinates(
				numberOfConformations,
				coordinatesPerConformation,
				allCoordinates);

		// Do the calculations
		calcRMSDOfOneVsFollowingWithRotation<<<this->blocks_per_grid, this->threads_per_block>>>(
				 this->deviceCoords,
				 reference_conformation_number,
				 reference_conformation_number+1,
				 numberOfConformations,
				 atomsPerConformation,
				 coordinatesPerConformation,
				 this->deviceRMSDs);
		
		// Update rotated coordinates
		updateHostCoordinates(
				numberOfConformations,
				coordinatesPerConformation,
				allCoordinates);
		
		if(rmsd!=NULL){
			// Get RMSDs
			updateHostRMSDs(
						numberOfConformations,
						reference_conformation_number,
						rmsd);
		}*/
}

void QCPCUDAKernel::oneVsFollowingFitDiffersCalcWithConfRotation(
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

