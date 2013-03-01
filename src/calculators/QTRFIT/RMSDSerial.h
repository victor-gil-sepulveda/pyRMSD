#ifndef RMSD_SERIAL_H_
#define RMSD_SERIAL_H_

#include <vector>
#include "../RMSD.h"

class RMSDSerial: public RMSD{
	
	public:
		RMSDSerial(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		~RMSDSerial();
	private:
		void _one_vs_following_fit_equals_calc_coords(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference, int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference, double* calcReference, int reference_conformation_number, double *rmsd);

};

#endif
