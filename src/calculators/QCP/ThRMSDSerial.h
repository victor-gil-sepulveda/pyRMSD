
#ifndef RMSD_CUDA_SERIAL_H_
#define RMSD_CUDA_SERIAL_H_

#include <vector>
#include "../RMSDCalculator.h"

class QCPSerialCalculator: public RMSDCalculator{

	public:
		QCPSerialCalculator(int numberOfConformations, int atomsPerConformation, double* coords);
		virtual ~QCPSerialCalculator();
		
	protected:

		KernelFunctions* getKernelFunctions();

};

#endif
