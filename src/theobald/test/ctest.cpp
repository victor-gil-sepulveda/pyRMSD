#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include "../ThRMSDSerial.h"
#include "../ThRMSDCuda.cuh"
using namespace std;

void load_vector(vector<double> & , const char * );
void print_vector(const char*,double*, int);
bool expectedVectorEqualsCalculatedWithinPrecision(const double * const , const double * const , int , double );
void compareVectors(const char*, const double * const , const double * const , int , double );

int main(int argc, char **argv){
	vector<double> serial_coordinates, cuda_coordinates, big_cuda_coordinates;
	vector<double> cuda_centered_coordinates;

	double* rmsd_serial = new double[11];
	double* rmsd_cuda = new double[11];
	vector<double> rmsd_cuda_big;

	// Load golden data
	load_vector(serial_coordinates, "data/coordinates");
	load_vector(cuda_coordinates, "data/coordinates");
	load_vector(big_cuda_coordinates, "data/coordinates_big");

	// Calculate
	ThRMSDSerial serial(11,66,&(serial_coordinates[0]));
	serial.oneVsTheOthers(0, rmsd_serial);
	print_vector("Serial:",rmsd_serial,11);

//	ThRMSDCuda cuda(11,66,&(cuda_coordinates[0]));
//	cuda.getDeviceCoordinates(cuda_centered_coordinates);
//	cuda.oneVsTheOthers(0, rmsd_cuda);
//	print_vector("CUDA:",rmsd_cuda,11);
//
//	// Assert!
//	compareVectors("Comparing centered coordinates (precision = 1e-4) ... ", serial.getCoordinates(),\
//			&(cuda_centered_coordinates[0]), 11*66*3, 1e-4);
//	compareVectors("Comparing rmsd_results (precision = 1e-4) ... ", &(rmsd_serial[0]),\
//			&(rmsd_cuda[0]), 11, 1e-4);

	// And clean the room...
	delete [] rmsd_serial;
	delete [] rmsd_cuda;


	// Now try with a big dataset (for profiling)
	ThRMSDCuda cuda_big(10001,66,&(big_cuda_coordinates[0]));
	cuda_big.calculateRMSDCondensedMatrix(rmsd_cuda_big);

	return 0;
}

void compareVectors(const char* message, const double * const expectedVector, const double * const calculatedVector, int dimension, double precision){
	cout<<message;
	bool comparison = expectedVectorEqualsCalculatedWithinPrecision(expectedVector,calculatedVector,dimension,precision);

	if(comparison == true){
		cout<<"OK"<<endl;
	}
	else{
		cout<<"KO"<<endl;
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

void print_vector(const char*message, double*  rmsd_vec, int len){
	cout<<message<<" ";
	for(int i =0; i< len; ++i){
		cout<<rmsd_vec[i]<<" ";
	}
	cout<<endl;
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
