
#ifndef RMSD_CUDA_SERIAL_H_
#define RMSD_CUDA_SERIAL_H_

#include <vector>
#include "../serial/RMSD.h"

#define floating_point_type double

class ThRMSDSerial: public RMSD{

	public:
		ThRMSDSerial(int numberOfConformations, int atomsPerConformation, double* coords);
		virtual ~ThRMSDSerial();
		virtual void oneVsFollowing(int conformation, double* rmsd_result);
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		
		// For debug purposes
		double* getCoordinates();

	protected:
		// Pointers to data storage
		floating_point_type* convertedCoords;
		floating_point_type* tmpRMSDs;
};

#endif
