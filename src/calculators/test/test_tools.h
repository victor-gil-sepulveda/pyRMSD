/*
 * test_tools.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef TEST_TOOLS_H_
#define TEST_TOOLS_H_
#include <vector>

void print_test_tittle(const char* function_name);
void load_vector(std::vector<double> & , const char * );
void print_vector(const char*,double*, int, int precission = 8);
bool expectedVectorEqualsCalculatedWithinPrecision(const double * const , const double * const , int , double );
void compareVectors(const char*, const double * const , const double * const , int , double );
void checkDistances(double* vector1, double* vector2, int totalatoms);
void writeVector(std::vector<double> & vector, const char* path);



#endif /* TEST_TOOLS_H_ */
