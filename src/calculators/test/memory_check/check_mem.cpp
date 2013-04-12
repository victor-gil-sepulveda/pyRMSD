#include "../test_tools.h"
#include "../../factory/RMSDCalculatorFactory.h"
#include "../../RMSDCalculator.h"

#include <iostream>
using namespace std;

void test_matrix_creation(RMSDCalculatorType type){
	print_test_tittle(__FUNCTION__);
	cout<<"- Using "<<calculatorTypeToString(type)<<endl;

	int number_of_coordsets = 5;
	int number_of_atoms = 224;

	vector<double> coordinates;
	vector<double> rmsds;

	load_vector(coordinates, "data/ligand_mini_CAs");

	test_vector_len(coordinates, number_of_atoms*number_of_coordsets*3, "CAs");

	// Iterposition
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
			type,
			number_of_coordsets,
			number_of_atoms,
			&(coordinates[0]));


	calculator->calculateRMSDCondensedMatrix(rmsds);

	delete calculator;
}

int main(int argc, char **argv){

	RMSDCalculatorType calculator_types [] =  {
			KABSCH_SERIAL_CALCULATOR,
			QTRFIT_SERIAL_CALCULATOR,
			QCP_SERIAL_CALCULATOR,
			KABSCH_OMP_CALCULATOR,
			QTRFIT_OMP_CALCULATOR,
			QCP_OMP_CALCULATOR
	};

	for(int i = 0; i< 6;++i) {
		test_matrix_creation(calculator_types[i]);
	}
}
