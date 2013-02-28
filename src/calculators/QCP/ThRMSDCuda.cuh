
#ifndef THEO_RMSD_CUDA_H_
#define THEO_RMSD_CUDA_H_

#include <vector>
#include "../RMSD.h"

#define floating_point_type float


class ThRMSDCuda: public RMSD{ 

	public:
		ThRMSDCuda(int numberOfConformations, int atomsPerConformation, double* coords, int threads_per_block = 32, int number_of_blocks = 512);
		~ThRMSDCuda();
		void oneVsFollowing(int conformation, double* rmsd_result );
		void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		
		// For debug purposes 
		void getDeviceCoordinates(std::vector<double>&);
		
	private:
		void cudaInit();
	
		// GPU-specific variables
		int numBlocks;
		int threadsPerBlock ;
		int numberOfBlocks;
		int rmsdMatrixLen;
	
		// Pointers to data storage
		floating_point_type* deviceCoords; 		// Pointer to device (GPU) data
		floating_point_type* deviceRMSDs;		// To store results
		floating_point_type* tmpHostRMSDs;    	// To store them on the host
		floating_point_type* tmpHostRMSDMatrix; // To store the RMSD Matrix
		floating_point_type* tmpHostCoords;    	// To store converted coordinates
};

#endif /* THEO_RMSD_CUDA_H_ */
