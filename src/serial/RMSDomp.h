/*
 * RMSDomp.h
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#ifndef RMSDOMP_H_
#define RMSDOMP_H_
#include "RMSD.h"


class RMSDomp: public RMSD{

	public:
		RMSDomp(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		virtual ~RMSDomp();
		void oneVsFollowing(int conformation, double* rmsd);
		void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		void iterativeSuperposition();
};

#endif /* RMSDOMP_H_ */
