
#ifndef RMSD_CUDA_H_
#define RMSD_CUDA_H_

#include <vector>
#include "../serial/RMSD.h"

class RMSDCuda: public RMSD{ 

	public:
		RMSDCuda(int numberOfConformations, int atomsPerConformation, double* coords, int threads_per_block = 64);
		~RMSDCuda();
		void oneVsTheOthers(int conformation, double* rmsd_result);
		void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		
	private:
	
		// GPU-specific variables
		int numBlocks;
		int threadsPerBlock ;
	
		// Pointers to data storage
		float* deviceCoords; 	// Pointer to device (GPU) data
		float* deviceRMSDs;		// To store results
		float* tmpHostRMSDs;    // To store them on the host
		float* coalescedCoords;	// Coordinates in coalesced form
};

#endif /* RMSD_H_ */
