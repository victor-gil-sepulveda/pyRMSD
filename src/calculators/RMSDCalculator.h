#ifndef RMSD_H_
#define RMSD_H_

#include <vector>
#include "KernelFunctions.h"
#include <cstddef>

struct RMSDCommandParameters{
	int numberOfConformations;
	int atomsPerConformation;
	double* allCoordinates;
};

/*
 * This is the base class for an RMSD Calculator, and holds almost all its logic (different calculators are
 * indeed different versions of the kernels).
 */
class RMSDCalculator{

	public:
		RMSDCalculator(int numberOfConformations, int atomsPerConformation, double* allCoordinates, KernelFunctions* kernelFunctions);
		virtual ~RMSDCalculator();

		void setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates);

		virtual void oneVsFollowing(int conformation, double* rmsd);
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd);
		void iterativeSuperposition(double rmsd_diff_to_stop = 1e-4, double* iteration_rmsd = NULL);
		double iterative_superposition_step(double* reference_coords, double* mean_coords);

		inline void setCoordinatesRotationTo(bool this_val){
			this->rotateFittingCoordinates =  this_val;
		}

	protected:
		void calculateRMSDCondensedMatrixWithFittingCoordinates(std::vector<double>& rmsd);

		void calculateRMSDCondensedMatrixWithFittingAndCalculationCoordinates(std::vector<double>& rmsd);

		void superposition_with_external_reference(double*);

		void superposition_with_external_reference_without_calc_coords(double*);

		void superposition_with_external_reference_rotating_calc_coords(double*);


		virtual void _one_vs_following_fit_equals_calc_coords(double* reference,
				int reference_conformation_number, double *rmsd);

		virtual void _one_vs_following_fit_differs_calc_coords(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd);

		virtual void _one_vs_following_fit_equals_calc_coords_rotating_coordinates(double* reference,
				int reference_conformation_number, double *rmsd);

		virtual void _one_vs_following_fit_differs_calc_coords_rotating_coordinates(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd);

		int numberOfConformations;
		int atomsPerFittingConformation;
		int coordinatesPerFittingConformation;
		double* allFittingCoordinates; // Coordinates for fitting and RMSD (if allRMSDCoordinates == NULL)
		int atomsPerCalculationConformation;
		int coordinatesPerCalculationConformation;
		double* allCalculationCoordinates; 	 // If is different from NULL, then this are the coordinates
									 	 	 // to calculate the RMSD.

		bool rotateFittingCoordinates;

		KernelFunctions* kernelFunctions;

	private:
		RMSDCalculator(){}
		RMSDCalculator(const RMSDCalculator& r){}
};
#endif
