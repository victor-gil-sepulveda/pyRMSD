/*
 * RMSDomp.h
 *
 *  Created on: 01/08/2012
 *      Author: victor
 */

#ifndef RMSDOMP_H_
#define RMSDOMP_H_
#include "../RMSDCalculator.h"


class QTRFITOmpCalculator: public RMSDCalculator{

	public:
		QTRFITOmpCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		virtual ~QTRFITOmpCalculator();

	private:
		void _one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);

		KernelFunctions* getKernelFunctions();
};

#endif /* RMSDOMP_H_ */
