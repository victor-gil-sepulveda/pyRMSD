#include "tests.h"
#include "../factory/RMSDCalculatorTypes.h"
#include <vector>
#include <cstdlib>
using std::vector;


void set_precision_if(
		RMSDCalculatorType calctype,
		double alternative_precision,
		double else_precision,
		double& precision){
	precision = else_precision;
	#ifdef USE_CUDA
		#ifdef CUDA_PRECISION_SINGLE
		if (calctype == QCP_CUDA_CALCULATOR)
			precision = alternative_precision;
		#endif
	#endif
}

int main(int argc, char **argv){

	RMSDCalculatorType available_calculators_d [] =  {
			KABSCH_SERIAL_CALCULATOR,
			KABSCH_OMP_CALCULATOR,
			QTRFIT_SERIAL_CALCULATOR,
			QTRFIT_OMP_CALCULATOR,
			QCP_SERIAL_CALCULATOR,
			QCP_OMP_CALCULATOR,
#ifdef USE_CUDA
			QCP_CUDA_CALCULATOR
#endif
	};

	vector<RMSDCalculatorType> available_calculators( available_calculators_d,
			available_calculators_d + sizeof(available_calculators_d)/sizeof(RMSDCalculatorType));

	test_initialize();
	test_copy_array();
	test_coordinates_mean();
	test_translations();

	test_QCP_Kernel();
	test_KABSCH_Kernel();

	// This tests does not work with QCP and KABSCH
	test_superposition_with_coordinates_change(QTRFIT_SERIAL_CALCULATOR);
	test_superposition_with_coordinates_change(QTRFIT_OMP_CALCULATOR);

	test_superposition_with_different_fit_and_calc_coordsets(QTRFIT_SERIAL_CALCULATOR, 1e-12);

	test_iterative_superposition_with_different_calc_and_fit_sets(QTRFIT_SERIAL_CALCULATOR, 1e-6);
	// Do those for all the others
/*	for(unsigned int i = 0; i < available_calculators.size();++i){
		double precision;

		test_superposition_with_very_different_calc_and_fit_sets(available_calculators[i], 1e-6);

		set_precision_if(
				available_calculators[i],
				1e-4,
				1e-6,
				precision);

		test_iterative_superposition_with_equal_calc_and_fit_sets(available_calculators[i], precision);

		set_precision_if(
					available_calculators[i],
					1e-4,
					1e-12,
					precision);
		test_superposition_with_different_fit_and_calc_coordsets(available_calculators[i], precision);

		set_precision_if(
					available_calculators[i],
					1e-4,
					1e-6,
					precision);
		test_iterative_superposition_with_different_calc_and_fit_sets(available_calculators[i], precision);
	}*/

	return 0;
}

