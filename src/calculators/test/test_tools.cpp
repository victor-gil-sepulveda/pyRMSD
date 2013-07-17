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
#include "../RMSDTools.h"
using namespace std;

void print_test_tittle(const char* function_name){
	cout <<"\n\033[34mTESTING "<<function_name<<"\033[0m"<<endl;
}

void print_calculator_and_precission(RMSDCalculatorType type, double precission){
	cout<<"- Using "<<calculatorTypeToString(type)<<" (prec. "<<setprecision(1)<<precission<<"):"<<endl;
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
	bool comparison = false;
	comparison = expectedVectorEqualsCalculatedWithinPrecision(expectedVector,calculatedVector,dimension,precision);

	if(comparison == true){
		cout<<"\033[1;32mOK\033[0m"<<endl;
	}
	else{
		cout<<"\033[1;31mKO\033[0m"<<endl;
	}
}

bool expectedVectorEqualsCalculatedWithinPrecision(
		const double * const expectedVector,
		const double * const calculatedVector,
		int dimension,
		double precision){
	bool equal = true;

    for(int i = 0; i < dimension; ++i){
    	bool is_nan = isnan(calculatedVector[i]);
    	bool has_big_error = fabs(expectedVector[i]-calculatedVector[i]) >= precision;
        if( is_nan || has_big_error ){
            equal = false;
            cout<<setprecision(16)<<"Problem: expectedVector["<<i<<"]="<<expectedVector[i]<<
            		" calculatedVector["<<i<<"]="<<calculatedVector[i]<<endl;
            cout<<" (dif = "<<fabs(expectedVector[i]-calculatedVector[i])<<")"<<endl;
            break;
        }
    }
    return equal;
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
		cout << "\033[1;31mUnable to open file: "<< string(dataPath) <<"\033[0m"<< endl;
	}
}

void load_pdb_coords(vector<double> & coords, vector<int> & shape, const char * dataPath){

	string line;
	ifstream myfile (dataPath);

	int numConfs, numAtoms, dim;
	double x,y,z;

	coords.clear();
	shape.clear();

	if (myfile.is_open()){
		getline(myfile, line);
		// First goes shape
		stringstream ss(line);
		ss >> numConfs >> numAtoms >> dim;
		shape.clear();
		shape.push_back(numConfs);
		shape.push_back(numAtoms);
		shape.push_back(dim);

		while(getline(myfile, line)){
			stringstream ss(line);
			ss >> x >> y >> z;
			coords.push_back(x);
			coords.push_back(y);
			coords.push_back(z);
		}

		myfile.close();
	}
	else{
		cout << "\033[1;31mUnable to open file: "<< string(dataPath) <<"\033[0m"<< endl;
	}
}

void load_and_center_pdb_coords(vector<double> & coords,
		vector<int> & shape,
		const char * dataPath,
		vector<double>* centers){

	load_pdb_coords(coords, shape, dataPath);

	if (centers!= NULL){
		centers->resize(shape[0]*3);
		RMSDTools::centerAllAtOrigin(
						shape[1],
						shape[0],
						&(coords[0]),
						&((*centers)[0])
		);
	}
	else{
		RMSDTools::centerAllAtOrigin(
								shape[1],
								shape[0],
								&(coords[0])
		);
	}
}

void load_and_move_pdb_coords(std::vector<double> & coords,
									std::vector<int> & shape,
									const char * dataPath,
									double* centers){
	load_pdb_coords(coords, shape, dataPath);
	RMSDTools::applyTranslationsToAll(
			(unsigned int) shape[1],
			(unsigned int) shape[0],
			&(coords[0]),
			centers,
			-1
	);
}

void save_pdb_coords(vector<double> & coords, vector<int> & shape, const char * file){
	ofstream ofs (file, std::ofstream::out);
	ofs.setf (std::ios::scientific);
	ofs<<shape[0]<<" "<<shape[1]<<" "<<shape[2]<<" "<<endl;
	for(unsigned int i = 0; i < coords.size(); i+=3){
		ofs<<coords[i]<<" "<<coords[i+1]<<" "<<coords[i+2]<<endl;
	}
	ofs.close();
}

void test_vector_len(vector<double>& v, unsigned int expected_len, const char* name){
	if(v.size()!= expected_len){
		cout<<"\033[1;31m"<<name<<" vector size is "<<v.size()<<" instead of "<<expected_len<<"\033[0m"<<endl;
	}
}

