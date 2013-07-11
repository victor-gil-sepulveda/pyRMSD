/*
 * tests.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef TESTS_H_
#define TESTS_H_
#include "../factory/RMSDCalculatorTypes.h"
	void test_initialize();
	void test_copy_array();
	void test_mean_coordinates();
	void test_translations();
	void test_center_coordinates();

	void test_QCP_Kernel();
	void test_KABSCH_Kernel();

	void test_superposition_with_fit(RMSDCalculatorType type,
			const char* initial_coords_file,
			const char* final_coords_file,
			const char* rmsd_results_file,
			double precision_of_check);

	void test_superposition_with_fit_and_calc(RMSDCalculatorType type,
			const char* initial_prot_coords_file,
			const char* final_prot_coords_file,
			const char* initial_lig_coords_file,
			const char* final_lig_coords_file,
			const char* rmsd_results_file,
			double precision_of_check);


#endif /* TESTS_H_ */
