/*
 * QCPCUDAMemKernel.cpp
 *
 *  Created on: Apr 13, 2013
 *      Author: victor
 */

#include "QCPCUDAMemKernel.h"
#include "kernel_functions_cuda.h"
inline void checkCudaError(char* message, cudaError error_code){
	if (error_code != 0){
		std::cout<<"Error in "<<message<<". Error code: "<<error_code<<". Exiting..."<<std::flush<<std::endl;
		exit(-1);
	}
}

QCPCUDAMemKernel::QCPCUDAMemKernel(
				double* coordinates,
				int atomsPerConformation,
				int coordinatesPerConformation,
				int numberOfConformations,
				int threads_per_block,
				int blocks_per_grid):QCPCUDAKernel( coordinates,
													atomsPerConformation,
													coordinatesPerConformation,
													numberOfConformations,
													threads_per_block,
													blocks_per_grid) {

}

QCPCUDAMemKernel::~QCPCUDAMemKernel(){}

void QCPCUDAMemKernel::matrixInit(
						double* allFittingCoordinates,
						int coordinatesPerFittingConformation,
						double* allCalculationCoordinates,
						int coordinatesPerCalculationConformation,
						int numberOfConformations){
	
	QCPCUDAKernel::matrixInit(
						allFittingCoordinates,
						coordinatesPerFittingConformation,
						allCalculationCoordinates,
						coordinatesPerCalculationConformation,
						numberOfConformations);
						
	//Allocate space to store all rmsds
	checkCudaError("Malloc allDeviceRMSDs", 
			cudaMalloc(
					(void **) &this->allDeviceRMSDs, 
					((numberOfConformations*(numberOfConformations-1)) / 2) * sizeof(floating_point_type)));
	
}

void QCPCUDAMemKernel::matrixEnd(double* rmsds_tmp, 
								int rmsds_tmp_len, 
								std::vector<double>& rmsds){
	
	
	
	rmsds.resize(rmsds_tmp_len);
	
	#ifdef CUDA_PRECISION_SINGLE
	
		float* buffer = new float[rmsds_tmp_len];
		
		checkCudaError("allDeviceRMSDs copy to host",
			cudaMemcpy(	buffer,
					this->allDeviceRMSDs,
					rmsds_tmp_len * sizeof(float),
					cudaMemcpyDeviceToHost));
		
		for(int i = 0; i < rmsds_tmp_len; ++i){
			rmsds[i] = static_cast<double>( buffer[i] );
		}
		
		delete [] buffer;
		
	#else	
		checkCudaError("allDeviceRMSDs copy to host",
			cudaMemcpy(	&(rmsds[0]),
					this->allDeviceRMSDs,
					rmsds_tmp_len * sizeof(double),
					cudaMemcpyDeviceToHost));
	#endif
	
	checkCudaError("allDeviceRMSDs cudaFree",
				cudaFree(this->allDeviceRMSDs));
}

int calculate_offset_for(int conformation_number, int number_of_conformations){
	/*int offset = 0;
	
	for(int i = 0; i < conformation_number; ++i){
		offset += (number_of_conformations-i-1);
	}
	
	return offset;*/
	return ((number_of_conformations-1)* conformation_number) - (conformation_number*(conformation_number-1))/2;
}

void QCPCUDAMemKernel::matrixOneVsFollowingFitEqualCalcWithoutConfRotation(
									double* reference, 
									int reference_conformation_number, 
									double* rmsd,
									int numberOfConformations, 
									int coordinatesPerConformation,
									int atomsPerConformation, 
									double* allCoordinates){
	
	floating_point_type* tmpDeviceReference = &(this->deviceCoords[reference_conformation_number*coordinatesPerConformation]);
	
	// Do the calculations
	calcRMSDOfOneVsFollowing CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			tmpDeviceReference,
			reference_conformation_number,
			this->deviceCoords,
			numberOfConformations,
			atomsPerConformation,
			coordinatesPerConformation,
			&(this->allDeviceRMSDs[calculate_offset_for(reference_conformation_number, numberOfConformations)]));
}

void QCPCUDAMemKernel::matrixOneVsFollowingFitDiffersCalcWithoutConfRotation(
											double* fitReference, 
											double* calcReference,
											int reference_conformation_number, 
											double* rmsd,
											int numberOfConformations, 
											int coordinatesPerConformation,
											int atomsPerConformation, 
											double* allCoordinates,
											int coordinatesPerRMSDConformation, 
											int atomsPerRMSDConformation,
											double* allRMSDCoordinates){
		
	floating_point_type* tmpFitDeviceReference = &(this->deviceCoords[reference_conformation_number*coordinatesPerConformation]);
	floating_point_type* tmpCalcDeviceReference = &(this->deviceCalcCoords[reference_conformation_number*coordinatesPerRMSDConformation]);
	
	// Do the calculations
	calcRMSDOfOneVsFollowingFitDiffersCalc CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			tmpFitDeviceReference,
			tmpCalcDeviceReference,
			reference_conformation_number,
			&(this->allDeviceRMSDs[calculate_offset_for(reference_conformation_number, numberOfConformations)]),
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			deviceCoords,
			coordinatesPerRMSDConformation,
			atomsPerRMSDConformation,
			deviceCalcCoords);
}
