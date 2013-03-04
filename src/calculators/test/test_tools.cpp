/*
 * test_tools.cpp
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#include "test_tools.h"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <iostream>
#include <cmath>
using namespace std;

void print_test_tittle(const char* function_name){
	cout <<"\n\033[34mTESTING "<<function_name<<"\033[0m"<<endl;
}

void checkDistances(double* vector1, double* vector2, int totalatoms){
	for(int i=0 ; i< totalatoms;i++){
		double* point1 = &(vector1[i*3]);
		double* point2 = &(vector2[i*3]);
		double tmp = 	(point1[0]-point2[0])*(point1[0]-point2[0]) +
						(point1[1]-point2[1])*(point1[1]-point2[1]) +
						(point1[2]-point2[2])*(point1[2]-point2[2]);
		cout<<sqrt(tmp)<<endl;
	}
}

void writeVector(vector<double> & vector, const char* path){
	std::ofstream output_file(path);
	output_file.precision(12);
	std::ostream_iterator<double> output_iterator(output_file, "\n");
	std::copy(vector.begin(), vector.end(), output_iterator);
}

void compareVectors(const char* message, const double * const expectedVector, const double * const calculatedVector, int dimension, double precision){
	cout<<message;
	bool comparison = expectedVectorEqualsCalculatedWithinPrecision(expectedVector,calculatedVector,dimension,precision);

	if(comparison == true){
		cout<<"\033[1;32mOK\033[0m"<<endl;
	}
	else{
		cout<<"\033[1;31mKO\033[0m"<<endl;
	}
}

bool expectedVectorEqualsCalculatedWithinPrecision(const double * const expectedVector, const double * const calculatedVector, int dimension, double precision){
	bool equal = true;

    for(int i=0;i<dimension;++i)
    {
        if( fabs(expectedVector[i]-calculatedVector[i]) >= precision )
        {
            equal = false;

            cout<<setprecision(16)<<"Problem: expectedVector["<<i<<"]="<<expectedVector[i]<<
            		" calculatedVector["<<i<<"]="<<calculatedVector[i]<<endl;
            cout<<" (dif = "<<fabs(expectedVector[i]-calculatedVector[i])<<")"<<endl;
            break;
        }
    }
    return equal;
}

void print_vector(const char*message, double*  rmsd_vec, int len, int precission){
	cout<<message<<" ";
	for(int i =0; i< len; ++i){
		cout<<setprecision(precission)<<rmsd_vec[i]<<" ";
	}
	cout<<flush<<endl;
}


inline double toDouble(const std::string & s){
	std::istringstream i(s);
	double x;
	i >> x;
	return x;
}

void load_vector(vector<double> & vector, const char * dataPath){
	string line;

	ifstream myfile (dataPath);

	if (myfile.is_open())	{
		while(getline(myfile, line)){
				vector.push_back(toDouble(line));
		}

		myfile.close();
	}
	else{
		cout << "Unable to open file: "<< string(dataPath) << endl;
	}
}



