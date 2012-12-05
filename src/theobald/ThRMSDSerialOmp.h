
#ifndef RMSD_CUDA_SERIAL_OMP_H_
#define RMSD_CUDA_SERIAL_OMP_H_

#include <vector>
#include "ThRMSDSerial.h"

#define floating_point_type double

class ThRMSDSerialOmp: public ThRMSDSerial{

	public:
		ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation, double* coords, int omp_threads);
		void oneVsFollowing(int conformation, double* rmsd_result);

	private:
		int omp_threads;
};

#endif
