#ifndef RMSD_SERIAL_H_
#define RMSD_SERIAL_H_

#include <vector>
#include "../RMSDCalculator.h"

class QTRFITSerialCalculator: public RMSDCalculator{
	
	public:
		QTRFITSerialCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		~QTRFITSerialCalculator();

	private:
		KernelFunctions* getKernelFunctions();

};

#endif
