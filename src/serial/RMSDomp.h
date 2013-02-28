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
		void superpositionChangingCoordinates(double*, double*);
		void iterativeSuperposition(double);

	private:
		void _one_vs_following_fit_equals_calc_coords(int conformation, double *rmsd, bool preserve_coords);
		void _one_vs_following_fit_differs_calc_coords(int conformation, double *rmsd, bool preserve_coords);
		void _one_vs_following_fit_equals_calc_coords_changing_coordinates(int conformation, double *rmsd, bool preserve_coords);
		void _one_vs_following_fit_differs_calc_coords_changing_coordinates(int conformation, double *rmsd, bool preserve_coords);
};

#endif /* RMSDOMP_H_ */
