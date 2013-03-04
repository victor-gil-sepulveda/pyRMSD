#include "tests.h"

int main(int argc, char **argv){

	test_initialize();
	test_copy_array();
	test_coordinates_mean();
	test_translations();
	test_iterative_superposition();
	test_iterative_superposition_with_different_calc_and_fit_sets();
	test_superposition_with_coordinates_change();
	test_superposition_with_different_fit_and_calc_coordsets();
	test_QCP_Kernel();

	return 0;
}

