#include <cmath>
#include <iostream>
#include "kernel_functions_cuda.h"
#include "QCPCUDAKernel.h"
#include "../RMSDTools.h"
using namespace std;

//
//	    0   0.0%   0.0%    31771  88.1% __libc_start_main
//      0   0.0%   0.0%    31199  86.5% calculateRMSDCondensedMatrix
//      3   0.0%   0.0%    31124  86.3% RMSDCalculator::calculateRMSDCondensedMatrix
//      1   0.0%   0.0%    30960  85.8% RMSDCalculator::_one_vs_following_fit_equals_calc_coords
//   8589  23.8%  23.8%    22851  63.3% cuMemGetAttribute_v2
//      4   0.0%  23.8%    21019  58.3% cudaGetExportTable
//      2   0.0%  23.8%    20993  58.2% cudaMemcpy
//      2   0.0%  23.8%    15926  44.2% QCPCUDAKernel::oneVsFollowingFitEqualCalcWithoutConfRotation
//   9864  27.3%  51.2%    15895  44.1% QCPCUDAKernel::updateDeviceCoordinates
//     73   0.2%  51.4%    15034  41.7% QCPCUDAKernel::updateHostRMSDs
//      0   0.0%  51.4%    14960  41.5% cuMemcpyDtoH_v2
//      1   0.0%  51.4%     6023  16.7% cuMemcpyHtoD_v2


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
inline void checkCudaError(char* message, cudaError error_code){
	if (error_code != 0){
		std::cout<<"Error in "<<message<<". Error code: "<<error_code<<". Exiting..."<<std::flush<<std::endl;
		exit(-1);
	}
}

void QCPCUDAKernel::init(
		double* coordinates, 
		int atomsPerConformation,
		int coordinatesPerConformation, 
		int numberOfConformations){

	
//	Now I'm  not using streams anymore as it didn't give me any speedup
//	int device;
//	cudaDeviceProp props;
//	
//	checkCudaError("cudaGetDevice",
//			cudaGetDevice(
//					&device));
//	
//	checkCudaError("cudaGetDeviceProperties",
//			cudaGetDeviceProperties(
//					&props, 
//					device));
//
//	if(!props.deviceOverlap){
//		cout<<"There are no multiple streams, so this implementation may not work properly :S"<<endl;
//	}

	int totalNumberOfCoordinates = coordinatesPerConformation*numberOfConformations;

	#ifdef CUDA_PRECISION_SINGLE
		// Allocate space for temporary coords
		this->tmpHostCoords = new float[totalNumberOfCoordinates];
		this->tmpHostRMSDs = new float[numberOfConformations];
		this->tmpHostReference = new float[coordinatesPerConformation];

	#endif
		
	// GPU Data allocation for input
	checkCudaError("Malloc Device Coords", 
				cudaMalloc(
						(void **) &deviceReference, 
						coordinatesPerConformation * sizeof(floating_point_type)));
	
	checkCudaError("Malloc Device Coords", 
			cudaMalloc(
					(void **) &deviceCoords, 
					totalNumberOfCoordinates * sizeof(floating_point_type)));
	

	// GPU Data allocation for output
	checkCudaError("Malloc RMSDs", 
			cudaMalloc(
					(void **) &deviceRMSDs, 
					numberOfConformations * sizeof(floating_point_type)));
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

	if (tmpHostReference != NULL){
		delete [] tmpHostReference;
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

	if (tmpCalcHostReference != NULL){
		delete [] tmpCalcHostReference;
	}

	if (deviceCalcReference != NULL){
		checkCudaError("deviceCalcReference cudaFree", 
				cudaFree(this->deviceCalcReference));
	}
}

void QCPCUDAKernel::changeCalculationCoords(
		double* calcCoords, 
		int number_of_atoms, 
		int numberOfConformations){
	
	#ifdef CUDA_PRECISION_SINGLE
			// Allocate space for temporary coords and copy contents
			this->tmpCalcHostReference = new float[number_of_atoms*3];
			this->tmpCalcHostCoords = new float[numberOfConformations*number_of_atoms*3];
	#endif
	
	checkCudaError("Malloc Device Calc Coords", 
				cudaMalloc(
						(void **) &deviceCalcReference, 
						number_of_atoms * 3 * sizeof(floating_point_type)));
	
	checkCudaError("Malloc Device Calc Coords", 
			cudaMalloc(
					(void **) &deviceCalcCoords, 
					number_of_atoms * 3 * numberOfConformations * sizeof(floating_point_type)));
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
	
	// Update reference in device
	updateDeviceCoordinates(
			reference,
			deviceReference,
			tmpHostCoords,
			coordinatesPerConformation,
			1);
	
	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			allCoordinates,
			deviceCoords,
			tmpHostCoords,
			coordinatesPerConformation,
			numberOfConformations);

	// Do the calculations
	calcRMSDOfOneVsFollowing CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			this->deviceReference, 
			reference_conformation_number,
			this->deviceCoords,
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
	
	// Update references in device
	updateDeviceCoordinates(
			fitReference,
			deviceReference,
			tmpHostCoords,
			coordinatesPerConformation,
			1);
	
	if(calcReference!=NULL){
		updateDeviceCoordinates(
				calcReference,
				deviceCalcReference,
				tmpCalcHostReference,
				coordinatesPerRMSDConformation,
				1);
	}
	
	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			allCoordinates,
			deviceCoords,
			tmpHostCoords,
			coordinatesPerConformation,
			numberOfConformations);

	updateDeviceCoordinates(
			allRMSDCoordinates,
			deviceCalcCoords,
			tmpCalcHostCoords,
			coordinatesPerRMSDConformation,
			numberOfConformations);

	// Do the calculations
	calcRMSDOfOneVsFollowingFitDiffersCalc CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			deviceReference,
			deviceCalcReference,
			reference_conformation_number,
			deviceRMSDs,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			deviceCoords,
			coordinatesPerRMSDConformation,
			atomsPerRMSDConformation,
			deviceCalcCoords);

	if(rmsd!=NULL){
		// Get RMSDs
		updateHostRMSDs(
					numberOfConformations,
					reference_conformation_number,
					rmsd);
				
	}
}

void QCPCUDAKernel::oneVsFollowingFitEqualCalcWithConfRotation(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		double *allCoordinates){
		
	// Update reference in device
	updateDeviceCoordinates(
			reference,
			deviceReference,
			tmpHostCoords,
			coordinatesPerConformation,
			1);
	
	// Put the centered coordinates on the device
	updateDeviceCoordinates(
			allCoordinates,
			deviceCoords,
			tmpHostCoords,
			coordinatesPerConformation,
			numberOfConformations);

	// Do the calculations
	calcRMSDOfOneVsFollowingWithRotation CUDA_KERNEL_DIM(this->blocks_per_grid, this->threads_per_block)(
			this->deviceReference, 
			reference_conformation_number,
			this->deviceCoords,
			numberOfConformations,
			atomsPerConformation,
			coordinatesPerConformation,
			this->deviceRMSDs);

	updateHostCoordinates(
			allCoordinates,
			deviceCoords,
			tmpHostCoords,
			numberOfConformations,
			coordinatesPerConformation);
	
	if(rmsd!=NULL){
		// Get RMSDs
		updateHostRMSDs(
					numberOfConformations,
					reference_conformation_number,
					rmsd);
	}
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

	oneVsFollowingFitDiffersCalcWithoutConfRotation(
			fitReference,
			calcReference,
			reference_conformation_number,
			rmsd,
			numberOfConformations,
			coordinatesPerConformation,
			atomsPerConformation,
			allCoordinates,
			coordinatesPerRMSDConformation,
			atomsPerRMSDConformation,
			allRMSDCoordinates);

	updateHostCoordinates(
				allCoordinates,
				deviceCoords,
				tmpHostCoords,
				numberOfConformations,
				coordinatesPerConformation);
	
	updateHostCoordinates(
				allRMSDCoordinates,
				deviceCalcCoords,
				tmpCalcHostCoords,
				numberOfConformations,
				coordinatesPerConformation);
	
}

