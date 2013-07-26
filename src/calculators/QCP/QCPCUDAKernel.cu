#include <cmath>
#include <iostream>
#include "kernel_functions_cuda.h"
#include "QCPCUDAKernel.h"
#include "../RMSDTools.h"
#include "../RMSDCalculationData.h"
using namespace std;

/**
 *	Convenience function for CUDA error handling. If it captures an error, writes a user-specified
 *  message, and exits the program.
 *
 *  \param 	message [In] Message to print if something went wrong (usually calling CUDA function name).
 *  \param 	error_code [In] CUDA error code (It must be 0 if everything is OK).
 */
inline void checkCudaError(char* message, cudaError error_code){
	if (error_code != 0){
		std::cout<<"Error in "<<message<<". Error code: "<<error_code<<". Exiting..."<<std::flush<<std::endl;
		exit(-1);
	}
}

/**
 * Kernel creator.
 */
QCPCUDAKernel::QCPCUDAKernel(
				RMSDCalculationData* data,
				int threads_per_block, 
				int blocks_per_grid){
	
	//	Note: I'm  not using streams anymore as it didn't give me any speedup
	this->threads_per_block = threads_per_block;
	this->blocks_per_grid = blocks_per_grid;
	
	tmpHostCoords = tmpHostRMSDs = tmpCalcHostCoords = NULL;
	
	deviceCoords =  deviceRMSDs = deviceReference =
	deviceCalcCoords = deviceCalcReference = NULL;

	#ifdef CUDA_PRECISION_SINGLE
		// Allocate space for buffers
		this->tmpHostCoords = new float[data->fittingCoordinatesLength];
		this->tmpHostRMSDs = new float[data->numberOfConformations];
	#endif
		
	// GPU Data allocation for input
	checkCudaError("Malloc Device Coords", 
				cudaMalloc(
						(void **) &deviceReference, 
						data->fittingConformationLength * sizeof(floating_point_type)));
	
	checkCudaError("Malloc Device Coords", 
			cudaMalloc(
					(void **) &deviceCoords, 
					data->fittingCoordinatesLength * sizeof(floating_point_type)));
	

	// GPU Data allocation for output
	checkCudaError("Malloc RMSDs", 
			cudaMalloc(
					(void **) &deviceRMSDs, 
					data->numberOfConformations * sizeof(floating_point_type)));
}

void QCPCUDAKernel::setCalculationCoords(RMSDCalculationData* data){

	#ifdef CUDA_PRECISION_SINGLE
			// Allocate space for temporary coords and copy contents
			this->tmpCalcHostCoords = new float[data->calculationCoordinatesLength];
	#endif

	checkCudaError("Malloc Device Calc Reference",
				cudaMalloc(
						(void **) &deviceCalcReference,
						data->calculationConformationLength * sizeof(floating_point_type)));

	checkCudaError("Malloc Device Calc Coords",
			cudaMalloc(
					(void **) &deviceCalcCoords,
					data->calculationCoordinatesLength * sizeof(floating_point_type)));
}

QCPCUDAKernel::~QCPCUDAKernel(){
	checkCudaError("deviceReference cudaFree",
				cudaFree(this->deviceReference));
	checkCudaError("deviceCoords cudaFree", 
			cudaFree(this->deviceCoords));
	checkCudaError("deviceRMSDs cudaFree", 
			cudaFree(this->deviceRMSDs));
	if (tmpHostCoords != NULL){
		delete [] tmpHostCoords;
	}
	if (tmpHostRMSDs != NULL){
		delete [] tmpHostRMSDs;
	}
	if (tmpCalcHostCoords != NULL){
		delete [] tmpCalcHostCoords;
	}
	if (deviceCalcCoords != NULL){
		checkCudaError("deviceCalcCoords cudaFree", 
				cudaFree(this->deviceCalcCoords));
	}
	if (deviceCalcReference != NULL){
		checkCudaError("deviceCalcReference cudaFree", 
				cudaFree(this->deviceCalcReference));
	}
}


void QCPCUDAKernel::updateDeviceCoordinates(
		double * coordinates,
		floating_point_type* device_coordinates,
		float* tmp_coordinates,
		int coordinates_per_conformation,
		int number_of_conformations){
	
	int total_number_of_coordinates = coordinates_per_conformation*number_of_conformations;
	
	#ifdef CUDA_PRECISION_SINGLE
		// Convert to our floating point type
		for(int i = 0; i < total_number_of_coordinates; ++i){
			tmp_coordinates[i] = static_cast<float>( coordinates[i] );
		}
		
		checkCudaError("Copying Coords to Device (single)", 
					cudaMemcpy(
							device_coordinates, 
							tmp_coordinates, 
							total_number_of_coordinates * sizeof(float),
							cudaMemcpyHostToDevice));
	#else	
		checkCudaError("Copying Coords to Device (double)", 
				cudaMemcpy(
						device_coordinates, 
						coordinates, 
						total_number_of_coordinates * sizeof(double),
						cudaMemcpyHostToDevice));
	#endif
}

void QCPCUDAKernel::updateHostCoordinates(
		double * coordinates,
		floating_point_type * device_coordinates,
		float * tmp_coordinates,
		int number_of_conformations,
		int coordinates_per_conformation
		){
	
	int total_number_of_coordinates = coordinates_per_conformation*number_of_conformations;
	
	#ifdef CUDA_PRECISION_SINGLE
		// Convert to our floating point type
		checkCudaError("Copying Coords to Host (single)", 
				cudaMemcpy(
						tmp_coordinates, 
						device_coordinates,
						total_number_of_coordinates * sizeof(float),
						cudaMemcpyDeviceToHost));
		
		// Conversion
		for(int i = 0; i < total_number_of_coordinates; ++i){
			coordinates[i] = static_cast<double>( tmp_coordinates[i]);
		}
		
	#else	
		checkCudaError("Copying Coords to Host (single)", 
				cudaMemcpy(
						coordinates, 
						device_coordinates, 
						total_number_of_coordinates * sizeof(double),
						cudaMemcpyDeviceToHost));
	#endif
}

void QCPCUDAKernel::updateHostRMSDs(
				int numberOfConformations,
				int reference_conformation_number,
				double* rmsd){
	
	#ifdef CUDA_PRECISION_SINGLE
		// Get RMSDs from device
		checkCudaError("Getting RMSDs from Device (single)", 
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
		checkCudaError("Getting RMSDs from Device (double)", 
				cudaMemcpy(
						rmsd,
						this->deviceRMSDs,
						(numberOfConformations - reference_conformation_number-1) * sizeof(double),
						cudaMemcpyDeviceToHost));
	#endif

}


void QCPCUDAKernel::oneVsFollowingFitEqualCalcCoords(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	// Update reference in device
	updateDeviceCoordinates(
			reference,
			deviceReference,
			tmpHostCoords,
			data->fittingConformationLength,
			1);
	
	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			data->fittingCoordinates,
			deviceCoords,
			tmpHostCoords,
			data->fittingConformationLength,
			data->numberOfConformations);

	// Do the calculations
	calcRMSDOfOneVsFollowingWithRotation CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			this->deviceReference, 
			reference_conformation_number,
			this->deviceCoords,
			data->numberOfConformations,
			data->atomsPerFittingConformation,
			data->fittingConformationLength,
			this->deviceRMSDs);

	updateHostCoordinates(
			data->fittingCoordinates,
			deviceCoords,
			tmpHostCoords,
			data->numberOfConformations,
			data->fittingConformationLength);
	
	if(rmsd!=NULL){
		// Get RMSDs
		updateHostRMSDs(
					data->numberOfConformations,
					reference_conformation_number,
					rmsd);
	}
}


void QCPCUDAKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	// Update references in device
	updateDeviceCoordinates(
			fitReference,
			deviceReference,
			tmpHostCoords,
			data->fittingConformationLength,
			1);

	// Update the (already centered) coordinates on the device
	updateDeviceCoordinates(
			data->fittingCoordinates,
			deviceCoords,
			tmpHostCoords,
			data->fittingConformationLength,
			data->numberOfConformations);

	if(calcReference!=NULL){
		updateDeviceCoordinates(
				calcReference,
				deviceCalcReference,
				tmpCalcHostCoords,
				data->calculationConformationLength,
				1);
	}

	updateDeviceCoordinates(
			data->calculationCoordinates,
			deviceCalcCoords,
			tmpCalcHostCoords,
			data->calculationConformationLength,
			data->numberOfConformations);

	// Do the calculations
	calcRMSDOfOneVsFollowingFitDiffersCalc CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			deviceReference,
			deviceCalcReference,
			reference_conformation_number,
			deviceRMSDs,
			data->numberOfConformations,
			data->fittingConformationLength,
			data->atomsPerFittingConformation,
			deviceCoords,
			data->calculationConformationLength,
			data->atomsPerCalculationConformation,
			deviceCalcCoords);

	updateHostCoordinates(
				data->fittingCoordinates,
				deviceCoords,
				tmpHostCoords,
				data->numberOfConformations,
				data->fittingConformationLength);
	
	updateHostCoordinates(
				data->calculationCoordinates,
				deviceCalcCoords,
				tmpCalcHostCoords,
				data->numberOfConformations,
				data->calculationConformationLength);
	
	if(rmsd!=NULL){
		// Get RMSDs
		updateHostRMSDs(
					data->numberOfConformations,
					reference_conformation_number,
					rmsd);
	}
}


void QCPCUDAKernel::matrixOneVsFollowingFitEqualCalc(
									double* reference, 
									int reference_conformation_number, 
									double* rmsd,
									RMSDCalculationData* data){
	
	floating_point_type* tmpDeviceReference = &(this->deviceCoords[reference_conformation_number*
	                                                               data->fittingConformationLength]);
	
	// Do the calculations
	calcRMSDOfOneVsFollowing CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			tmpDeviceReference, 
			reference_conformation_number,
			this->deviceCoords,
			data->numberOfConformations,
			data->atomsPerFittingConformation,
			data->fittingConformationLength,
			this->deviceRMSDs);

	// Get RMSDs
	updateHostRMSDs(
			data->numberOfConformations,
			reference_conformation_number,
			rmsd);
}

void QCPCUDAKernel::matrixOneVsFollowingFitDiffersCalc(
											double* fitReference, 
											double* calcReference,
											int reference_conformation_number, 
											double* rmsd,
											RMSDCalculationData* data){
		
	floating_point_type* tmpFitDeviceReference = &(this->deviceCoords[reference_conformation_number*
	                                                                  data->fittingConformationLength]);

	floating_point_type* tmpCalcDeviceReference = &(this->deviceCalcCoords[reference_conformation_number*
	                                                                       data->calculationConformationLength]);
	
	// Do the calculations
	calcRMSDOfOneVsFollowingFitDiffersCalc CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			tmpFitDeviceReference,
			tmpCalcDeviceReference,
			reference_conformation_number,
			deviceRMSDs,
			data->numberOfConformations,
			data->fittingConformationLength,
			data->atomsPerFittingConformation,
			deviceCoords,
			data->calculationConformationLength,
			data->atomsPerCalculationConformation,
			deviceCalcCoords);

	updateHostRMSDs(
				data->numberOfConformations,
				reference_conformation_number,
				rmsd);
					
}

void QCPCUDAKernel::matrixInit(RMSDCalculationData* data){

	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			data->fittingCoordinates,
			deviceCoords,
			tmpHostCoords,
			data->fittingConformationLength,
			data->numberOfConformations);

	// And if we have the others, then update them
	if(data->hasCalculationCoordinatesSet()){
		updateDeviceCoordinates(
				data->calculationCoordinates,
				deviceCalcCoords,
				tmpCalcHostCoords,
				data->calculationConformationLength,
				data->numberOfConformations);
	}
}
