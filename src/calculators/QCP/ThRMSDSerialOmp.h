
#ifndef RMSD_CUDA_SERIAL_OMP_H_
#define RMSD_CUDA_SERIAL_OMP_H_

#include <vector>
#include "ThRMSDSerial.h"

class ThRMSDSerialOmp: public ThRMSDSerial{

	public:
		ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation, double* coords, int omp_threads);
		~ThRMSDSerialOmp();

	private:
		int omp_threads;
};

#endif
