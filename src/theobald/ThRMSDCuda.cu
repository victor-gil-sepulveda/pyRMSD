#include "ThRMSDCuda.cuh"
#include "kernel_functions_cuda.cuh"
#include <iostream>
using namespace std;

#define floating_point_type float

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
		cout<<"Error in "<<message<<" . Error code: "<<error_code<<". Exiting..."<<flush<<endl;
		exit(-1);
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
///	This function encapsulates initial CUDA queries to the GPU.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDCuda::cudaInit(){
	cudaDeviceProp props;
	int device;
	checkCudaError("cudaGetDevice",cudaGetDevice(&device));
	checkCudaError("cudaGetDeviceProperties",cudaGetDeviceProperties(&props, device));
	
	if(!props.deviceOverlap){
		cout<<"There are no multiple streams, so this implementation may not work propperly :S"<<endl;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class constructor. Allocates memory both in host and device, and initializes them.
///
/// \param 	numberOfConformations [In] Total number of conformations in the coordinates array.
///
/// \param 	atomsPerConformation [In] Number of atoms of every conformation.
///
/// \param 	coords [In] Coordinates array with numberOfConformations*atomsPerConformation*3 elements.
///
/// \param 	threads_per_block [In] Threads per block to be used in CUDA kernel calls.
///
/// \param 	number_of_blocks [In] Number of blocks to be used in CUDA kernel calls.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
ThRMSDCuda::ThRMSDCuda(int numberOfConformations, int atomsPerConformation, double* coords, int threads_per_block, int number_of_blocks):
				RMSD(numberOfConformations, atomsPerConformation, coords){
    
    cudaInit();
    
	int total_num_of_coords = this->numberOfConformations*this->atomsPerConformation*3;
    
    // Convert to our floating point type
    this->tmpHostCoords = new floating_point_type[total_num_of_coords];
    for(int i = 0; i < total_num_of_coords; ++i){
    	this->tmpHostCoords[i] = (floating_point_type) coords[i];
    }
    
	// Allocate space for temporary rmsds
	this->tmpHostRMSDs = new floating_point_type[numberOfConformations];
    
	//////////////////////////////////////////////////////////
	//// Set up GPU
	//////////////////////////////////////////////////////////
	this->threadsPerBlock = threads_per_block;
	this->numberOfBlocks = number_of_blocks;
	
	
	// GPU Data allocation for input
	checkCudaError(" Malloc Device Coords. ", cudaMalloc((void **) &deviceCoords, total_num_of_coords * sizeof(floating_point_type)));
	checkCudaError(" Copying Coords to Device. ", cudaMemcpy(deviceCoords, tmpHostCoords, total_num_of_coords * sizeof(floating_point_type), 
			cudaMemcpyHostToDevice));
    
    // GPU Data allocation for output
    checkCudaError(" Malloc RMSDs ", cudaMalloc((void **) &deviceRMSDs, numberOfConformations * sizeof(floating_point_type)));
   
    
    
    //////////////////////////////////////////////////////////
	//// And do some pre-processing
	//////////////////////////////////////////////////////////
    centerCoordsOfAllConformations<<<numberOfBlocks, threadsPerBlock>>>(numberOfConformations,atomsPerConformation,deviceCoords);
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Class destructor. Frees memory.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
ThRMSDCuda::~ThRMSDCuda(){
    checkCudaError(" deviceCoords cudaFree ", cudaFree(this->deviceCoords));
    checkCudaError(" deviceRMSDs cudaFree ", cudaFree(this->deviceRMSDs));
    delete [] this->tmpHostRMSDs;
    delete [] this->tmpHostCoords;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Implements a 'i' vs [i+1,M] rmsd conformation calculation, where M is the number of conformations.
///
/// \param 	conformation [In] Total number of conformations in the coordinates array.
///
/// \param 	rmsd_result [In/Out] Array where the rmsds will be stored (it needs to have the correct size!).
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDCuda::oneVsFollowing(int conformation, double* rmsd_result) {
	if (conformation < numberOfConformations){
	    calcRMSDOfOneVsOthers<<<numberOfBlocks,threadsPerBlock>>>(this->deviceCoords, conformation, conformation + 1, 
	    															numberOfConformations, atomsPerConformation, 
	    															atomsPerConformation*3, this->deviceRMSDs);
	    
	    checkCudaError(" Getting RMSDs from Device. ", cudaMemcpy(this->tmpHostRMSDs, this->deviceRMSDs, 
	    											   numberOfConformations * sizeof(floating_point_type), 
	    											   cudaMemcpyDeviceToHost));
	    
	    // Do the copy to the output vector (needs to have the correct size)
		int j = 0;
		for (int i = conformation + 1; i < numberOfConformations;++i,++j){
			rmsd_result[j] = (double) this->tmpHostRMSDs[i];
		}
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Implements the calculation of the upper triangular matrix in row major format of the rmsd distance between all
/// conformations stored in the coordinates.
///
/// \param 	rmsd [In/Out] Vector where the resulting RMSDs will be stored.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDCuda::calculateRMSDCondensedMatrix(vector<double>& rmsd){ 
	
	cudaStream_t copy_stream, execution_stream;
	cudaStreamCreate(&copy_stream);
	cudaStreamCreate(&execution_stream);
	
	this->rmsdMatrixLen = (numberOfConformations*(numberOfConformations-1))/2;
	//checkCudaError("tmpHostRMSD cudaHostAlloc ", cudaHostAlloc((void**)&(this->tmpHostRMSDMatrix), 
	//											rmsdMatrixLen*sizeof(floating_point_type),cudaHostAllocDefault));
	this->tmpHostRMSDMatrix = new floating_point_type[rmsdMatrixLen];
	
	float time;
	cudaEvent_t start;
	cudaEventCreate(&start);
	cudaEvent_t end;
	cudaEventCreate( &end);

	
		cudaEventRecord(start, 0);
		int rmsdMatrixOffset = 0;
		int numberOfCalculatedRmsds = 0;
		for(int conformation_number = 0;conformation_number < numberOfConformations; ++conformation_number){ 
			numberOfCalculatedRmsds = numberOfConformations-(conformation_number+1);
	    	
	    	calcRMSDOfOneVsOthers<<<numberOfBlocks,
	    						threadsPerBlock,
	    						0,execution_stream>>>(deviceCoords, conformation_number, conformation_number + 1, 
	    																					numberOfConformations, atomsPerConformation, 
	    																					atomsPerConformation*3, deviceRMSDs);
			cudaStreamSynchronize(execution_stream);
		
	    	checkCudaError(" Getting RMSDs from Device. ", cudaMemcpyAsync(	&(tmpHostRMSDMatrix[rmsdMatrixOffset]), 
	    																&(deviceRMSDs[conformation_number+1]), 
	    											   					numberOfCalculatedRmsds * sizeof(floating_point_type), 
	    											   					cudaMemcpyDeviceToHost,
	    											   					copy_stream));
	    	cudaStreamSynchronize(copy_stream);
	    	
			rmsdMatrixOffset += numberOfCalculatedRmsds;
		}
		
		for (int i = 0; i < rmsdMatrixLen; ++i){
			rmsd.push_back((double)this->tmpHostRMSDMatrix[i]);
		}
		
		cudaEventRecord(end, 0 );
	    cudaEventSynchronize(end);
	    
	    cudaEventElapsedTime(&time, start, end);	
	    cout<<"Time for calculations (ms):"<<time<<endl;	
	     
    cudaEventDestroy(start);
    cudaEventDestroy(end);
    
    cudaStreamDestroy(copy_stream);
	cudaStreamDestroy(execution_stream);
    //checkCudaError("tmpHostRMSD cudaFreeHost ",cudaFreeHost(tmpHostRMSDMatrix));
    delete [] this->tmpHostRMSDMatrix;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	Returns the coordinates stored in the device (for testing purposes).
///
/// \param 	coordinates [In/Out] Vector where the device coordinates will be .
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void ThRMSDCuda::getDeviceCoordinates(vector<double>& coordinates){
	int number_of_coords = numberOfConformations*atomsPerConformation*3;
    floating_point_type* tmpCoords = new floating_point_type[number_of_coords];
	checkCudaError(" Getting Coordinates from Device. ", 
		cudaMemcpy(tmpCoords, deviceCoords, number_of_coords * sizeof(floating_point_type), cudaMemcpyDeviceToHost)
	);
	coordinates.clear();
	for (int i = 0; i < number_of_coords; ++i){
		coordinates.push_back((double)tmpCoords[i]);
	}
	delete [] tmpCoords;
}