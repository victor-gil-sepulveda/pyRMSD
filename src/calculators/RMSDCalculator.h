#ifndef RMSD_H_
#define RMSD_H_

#include <vector>
#include "KernelFunctions.h"
#include <cstddef>

class RMSDCalculator{

	public:
		RMSDCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates, KernelFunctions* kernelFunctions);
		virtual ~RMSDCalculator();

		void setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates);

		virtual void oneVsFollowing(int conformation, double* rmsd);
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		void iterativeSuperposition(double rmsd_diff_to_stop = 1e-4);

		void superposition_with_external_reference_and_fit_equals_calc(double*, double*);
		void superposition_with_external_reference_and_fit_differs_calc(double*);

		inline void setCoordinatesModificationTo(bool this_val){
			this->modifyFittingCoordinates =  this_val;
		}

	protected:
		virtual void _one_vs_following_fit_equals_calc_coords(double* reference,
				int reference_conformation_number, double *rmsd, bool center = true);

		virtual void _one_vs_following_fit_differs_calc_coords(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd, bool center = true);

		virtual void _one_vs_following_fit_equals_calc_coords_changing_coordinates(double* reference,
				int reference_conformation_number, double *rmsd);

		virtual void _one_vs_following_fit_differs_calc_coords_changing_coordinates(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd);

		// Coordinates info
		int numberOfConformations;

		int atomsPerConformation;
		int coordinatesPerConformation;
		double* allCoordinates; // Coordinates for fitting and RMSD (if allRMSDCoordinates == NULL)

		int atomsPerRMSDConformation;
		int coordinatesPerRMSDConformation;
		double* allRMSDCoordinates; // If is different from NULL, then this are the coordinates
									 // to calculate the RMSD.

		bool modifyFittingCoordinates;

		KernelFunctions* kernelFunctions;

	private:
		RMSDCalculator(){}
		RMSDCalculator(const RMSDCalculator& r){}
};
#endif
