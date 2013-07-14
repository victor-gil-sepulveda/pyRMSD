/*
 * test_tools.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef TEST_TOOLS_H_
#define TEST_TOOLS_H_
#include <vector>
#include <stdlib.h>
#include <sstream>
#include "../factory/RMSDCalculatorTypes.h"
#include <iomanip>
#include <iostream>

using std::setprecision;

void print_test_tittle(const char* function_name);
void print_calculator_and_precission(RMSDCalculatorType type, double precission);
void load_vector(std::vector<double> & , const char * );
void print_vector(const char*,double*, int, int precission = 8);
bool expectedVectorEqualsCalculatedWithinPrecision(const double * const , const double * const , int , double );
void compareVectors(const char*, const double * const , const double * const , int , double );
void checkDistances(double* vector1, double* vector2, int totalatoms);
void writeVector(std::vector<double> & vector, const char* path);
void test_vector_len(std::vector<double>& v, unsigned int expected_len, const char* name);

void load_pdb_coords(std::vector<double> & coords, std::vector<int> & shape, const char * dataPath);
void load_and_center_pdb_coords(std::vector<double> & coords,
									std::vector<int> & shape,
									const char * dataPath,
									std::vector<double>* centers = NULL);
void load_and_move_pdb_coords(std::vector<double> & coords,
									std::vector<int> & shape,
									const char * dataPath,
									double* centers);

void save_pdb_coords(std::vector<double> & coords, std::vector<int> & shape, const char * file);

template <class T>
inline std::string toString(T data){
	return static_cast<std::ostringstream*>( &(std::ostringstream() << data) )->str();
}

template <class T>
void print_vector(const char*message, T*  rmsd_vec, int len, int precission){
	std::cout<<message<<" ";
	for(int i =0; i< len; ++i){
		std::cout<<setprecision(precission)<<rmsd_vec[i]<<" ";
	}
	std::cout<<std::flush<<std::endl;
}


#endif /* TEST_TOOLS_H_ */
