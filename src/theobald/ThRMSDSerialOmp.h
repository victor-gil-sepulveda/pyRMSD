
#ifndef RMSD_CUDA_SERIAL_OMP_H_
#define RMSD_CUDA_SERIAL_OMP_H_

#include <vector>
#include "ThRMSDSerial.h"

#define floating_point_type double

class ThRMSDSerialOmp: public ThRMSDSerial{

	public:
		ThRMSDSerialOmp(int numberOfConformations, int atomsPerConformation, double* coords);
		void oneVsFollowing(int conformation, double* rmsd_result);
};

#endif
