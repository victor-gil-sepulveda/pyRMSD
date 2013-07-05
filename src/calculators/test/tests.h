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
	void test_coordinates_mean();
	void test_translations();
	void test_superposition_with_very_different_calc_and_fit_sets(RMSDCalculatorType type, double precision_of_check);
	void test_iterative_superposition_with_equal_calc_and_fit_sets(RMSDCalculatorType type, double precision_of_check);
	void test_iterative_superposition_with_different_calc_and_fit_sets(RMSDCalculatorType type, double precision_check);
	void test_superposition_with_coordinates_change(RMSDCalculatorType type);
	void test_superposition_with_different_fit_and_calc_coordsets(RMSDCalculatorType type, double precision_check);
	void test_QCP_Kernel();
	void test_KABSCH_Kernel();
#endif /* TESTS_H_ */
