#ifndef RMSD_SERIAL_H_
#define RMSD_SERIAL_H_

#include <vector>
#include "RMSD.h"

class RMSDSerial: public RMSD{
	
	public:
		RMSDSerial(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		~RMSDSerial();
		void oneVsTheOthers(int conformation, double* rmsd);
		void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
};

#endif
