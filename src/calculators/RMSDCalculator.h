#ifndef RMSD_H_
#define RMSD_H_

#include <vector>
#include "KernelFunctions.h"
#include <cstddef>

class RMSDCalculationData;

/**
 * This class implements an RMSD Calculator, and holds almost all the logic to perform any
 * of the RMSD operations. The associated calculators kernel defined the way the superposition
 * is done as well as the level of parallelization.
 */
class RMSDCalculator{

	public:
						RMSDCalculator(RMSDCalculationData* rmsdData, KernelFunctions* kernelFunctions);
		virtual 		~RMSDCalculator();

		virtual void 	oneVsFollowing(int conformation, double* rmsd);
		virtual void 	calculateRMSDCondensedMatrix(std::vector<double>& rmsd);

		void 			iterativeSuperposition(double rmsd_diff_to_stop = 1e-4, double* iteration_rmsd = NULL);
		double 			iterative_superposition_step(double* reference_coords, double* mean_coords);

	protected:

		void calculate_rmsd_condensed_matrix_with_fitting_coordinates(std::vector<double>& rmsd);
		void calculate_rmsd_condensed_matrix_with_fitting_and_calculation_coordinates(std::vector<double>& rmsd);

		void superposition_with_external_reference(double*);
		void superposition_with_external_reference_without_calc_coords(double*);
		void superposition_with_external_reference_rotating_calc_coords(double*);

		void _one_vs_following_fit_equals_calc_coords(double* reference,
				int reference_conformation_number, double *rmsd);
		void _one_vs_following_fit_differs_calc_coords(double* fitReference,
				double* calcReference, int reference_conformation_number, double *rmsd);


		RMSDCalculationData* rmsdData;

		KernelFunctions* kernelFunctions;

	private:
		RMSDCalculator(){}
		RMSDCalculator(const RMSDCalculator& r){}
};
#endif
