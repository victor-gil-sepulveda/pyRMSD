
#ifndef RMSD_CUDA_SERIAL_OMP_H_
#define RMSD_CUDA_SERIAL_OMP_H_

#include <vector>
#include "ThRMSDSerial.h"

class QCPOmpCalculator: public QCPSerialCalculator{

	public:
		QCPOmpCalculator(int numberOfConformations, int atomsPerConformation, double* coords, int omp_threads);
		~QCPOmpCalculator();

	private:
		int omp_threads;
};

#endif
