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
		void oneVsFollowing(int reference_conformation_number, double* rmsd);
		void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		void superposition_with_external_reference(double*, double*);
		void iterativeSuperposition(double rmsd_diff_to_stop = 1e-4);

	private:
		void _one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);
};

#endif /* RMSDOMP_H_ */
