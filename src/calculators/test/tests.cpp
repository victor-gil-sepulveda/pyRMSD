#include "tests.h"
#include "test_tools.h"
#include <vector>
#include <iostream>
#include "../RMSDTools.h"
#include "../RMSDCalculator.h"
#include "../QCP/QCPSerialKernel.h"
#include "../QTRFIT/QTRFITOmpKernel.h"
#include "../factory/RMSDCalculatorFactory.h"
#include "../KABSCH/KABSCHSerialKernel.h"

using namespace std;

void test_initialize(){
	print_test_tittle(__FUNCTION__);
	double expected_initialized[] = {3,3,3,3,3,3,3,3,3,3};
	double initialiazed[10];
	double matrix_mode[3][3];
	double matrix_mode_copy[9];
	double expected_matrix_mode_copy[] ={5,5,5,5,5,5,5,5,5};

	RMSDTools::initializeTo(initialiazed, 3, 10);
	compareVectors("\tTesting initialization: ", expected_initialized, initialiazed, 10, 1e-16);

	RMSDTools::initializeTo(matrix_mode[0],5, 9);
	int k = 0;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			matrix_mode_copy[k]=matrix_mode[i][j];
			k = k+1;
		}
	}
	compareVectors("\tTesting initialization (matrix): ", expected_matrix_mode_copy, matrix_mode_copy, 5, 1e-16);
}

void test_copy_array(){
	print_test_tittle(__FUNCTION__);
	double expected_array[] = {0,1,2,3,4,5,1234,7,8,9};
	double uninitialized_array[10];
	RMSDTools::copyArrays(uninitialized_array, expected_array, 10);
	compareVectors("\tTesting array copy: ", expected_array, uninitialized_array, 10, 1e-10);
}

void test_coordinates_mean(){
	print_test_tittle(__FUNCTION__);
	double mean_coordinates[9];
	double coordinates[] = { 1,1,1,  2,2,2,  3,3,3,  // 1st conformation
							  4,4,4,  5,5,5,  6,6,6,  // 2nd conformation
							  7,7,7,  8,8,8,  9,9,9 };// 3rd conformation
	double expected_mean_coordinates[] = {4, 4, 4, 5, 5, 5, 6, 6, 6};
	int number_of_atoms = 3;
	int number_of_conformations = 3;
	RMSDTools::calculateMeanCoordinates(mean_coordinates, coordinates,
				number_of_conformations, number_of_atoms);

	compareVectors("\tTesting coordinates mean: ", expected_mean_coordinates, mean_coordinates, 9, 1e-10);
}

void test_translations(){
	print_test_tittle(__FUNCTION__);
	double coordinates[] = {
			1,1,0,  2,2,0,  3,3,0,  // 1st conformation
			6,6,0,  5,5,0,  4,4,0,  // 2nd conformation (upside down)
			7,7,0,  8,8,0,  9,9,0,  // 3rd conformation
			5,4,0,  5,5,0,  4,4,0   // 4th conformation
	};

	double translations[] = {
			0,0,1,
			0,1,0,
			1,0,0,
			1,0,1
	};

	double expected_coordinates[]={
			1,1,1,  2,2,1,  3,3,1,
			6,7,0,  5,6,0,  4,5,0,
			8,7,0,  9,8,0,  10,9,0,
			6,4,1,  6,5,1,  5,4,1
	};

	double translation_vector [] = {1, 2, 4};

	double expected_retranslated_coordinates[]={
				2,3,5,  3,4,5,    4,5,5,
				7,9,4,  6,8,4,    5,7,4,
				9,9,4,  10,10,4,  11,11,4,
				7,6,5,  7,7,5,    6,6,5
	};

	RMSDTools::applyTranslationsToAll(3,4,coordinates,translations);

	compareVectors("\tTesting translated coordinates: ", expected_coordinates, coordinates, 3*3*4, 1e-12);

	RMSDTools::applyTranslationToAll(3,4,coordinates,translation_vector);

	compareVectors("\tTesting translated coordinates: ", expected_retranslated_coordinates, coordinates, 3*3*4, 1e-12);
}

void test_superposition_with_coordinates_change(RMSDCalculatorType type){
	print_test_tittle(__FUNCTION__);
	cout<<"- Using "<<calculatorTypeToString(type)<<":"<<endl;

	double coordinates[] = { 1,1,0,  2,2,0,  3,3,0,  // 1st conformation
							  6,6,0,  5,5,0,  4,4,0,  // 2nd conformation (upside down)
							  7,7,0,  8,8,0,  9,9,0,  // 3rd conformation
							  5,4,0,  5,5,0,  4,4,0}; // 4th conformation (upside down + one changed point
													  // with dist d((6,6),(5,4))
	double reference_coordinates[] = {4,4,0,  5,5,0,  6,6,0};
	double rmsds[] = {0,0,0,0};
	int number_of_atoms = 3;
	int number_of_conformations = 4;

	double expected_reference[] = {4,4,0,  5,5,0,  6,6,0};
	double expected_rmsds[] = {0, 0, 0, 0.91376624};
	double expected_coordinates[] = {4,4,0,  5,5,0,  6,6,0,
									  4,4,0,  5,5,0,  6,6,0,
			                          4,4,0,  5,5,0,  6,6,0,
			                          4.5286,5,0,  5.2357,4.29289,0,  5.2357,5.70711,0};

	double coordinates_copy [number_of_atoms*number_of_conformations*3];

	RMSDTools::copyArrays(coordinates_copy,coordinates,number_of_atoms*number_of_conformations*3);
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
			type,
			number_of_conformations,
			number_of_atoms,
			coordinates_copy);

	calculator->superposition_with_external_reference_and_fit_equals_calc(reference_coordinates, rmsds);

	compareVectors("\tTesting RMSD: ", expected_rmsds, rmsds, number_of_conformations, 1e-8);
	compareVectors("\tTesting coordinates: ", expected_coordinates, coordinates_copy,
				number_of_atoms*3*number_of_conformations, 1e-5);
	compareVectors("\tTesting reference: ", expected_reference, reference_coordinates, number_of_atoms*3, 1e-10);
	delete calculator;
}

void test_iterative_superposition_with_equal_calc_and_fit_sets(RMSDCalculatorType type, double precision_of_check){
	print_test_tittle(__FUNCTION__);
	cout<<"- Using "<<calculatorTypeToString(type)<<" (prec. "<<precision_of_check<<"):"<<endl;

	int number_of_coordsets = 5;
	int number_of_atoms = 3239;

	vector<double> not_aligned_coordinates;
	vector<double> iterposed_coordinates;

	load_vector(not_aligned_coordinates, "data/ligand_mini_all");
	load_vector(iterposed_coordinates, "data/ligand_mini_iterposed_all");

	test_vector_len(not_aligned_coordinates, number_of_atoms*number_of_coordsets*3, "all not aligned atoms");
	test_vector_len(iterposed_coordinates, number_of_atoms*number_of_coordsets*3, "all iterposed atoms");

	// Iterposition with QTRFIT
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
			type,
			number_of_coordsets,
			number_of_atoms,
			&(not_aligned_coordinates[0]));


	calculator->iterativeSuperposition();
	compareVectors("\tTesting iterposed coordinates: ",
			&(iterposed_coordinates[0]),
			&(not_aligned_coordinates[0]),
			not_aligned_coordinates.size(),
			precision_of_check);
	delete calculator;
}

void test_iterative_superposition_with_different_calc_and_fit_sets(RMSDCalculatorType type, double precision_check){
	print_test_tittle(__FUNCTION__);
	cout<<"- Using "<<calculatorTypeToString(type)<<" (prec. "<<precision_check<<"):"<<endl;

	int number_of_coordsets = 5;
	int number_of_atoms = 3239;
	int number_of_CAs = 224;

	vector<double> not_iterposed_coordinates;
	vector<double> not_aligned_CA;
	vector<double> iterposed_coordinates;

	load_vector(not_aligned_CA, "data/ligand_mini_CAs");
	load_vector(not_iterposed_coordinates, "data/ligand_mini_all");
	load_vector(iterposed_coordinates, "data/ligand_mini_iterposed_with_cas_all_atom");

	test_vector_len(not_aligned_CA, number_of_CAs*number_of_coordsets*3, "not aligned CAs");
	test_vector_len(not_iterposed_coordinates, number_of_atoms*number_of_coordsets*3, "all not aligned atoms");
	test_vector_len(iterposed_coordinates, number_of_atoms*number_of_coordsets*3, "all iterposed atoms");

	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
				type,
				number_of_coordsets,
				number_of_CAs,
				&(not_aligned_CA[0]));

	calculator->setCalculationCoordinates(number_of_atoms, &(not_iterposed_coordinates[0]));
	calculator->iterativeSuperposition();
	compareVectors("\tTesting iterposed coordinates: ", &(iterposed_coordinates[0]),&(not_iterposed_coordinates[0]), not_iterposed_coordinates.size(), precision_check);
	//writeVector( not_iterposed_coordinates, "vector.out");
	delete calculator;
}

void test_superposition_with_different_fit_and_calc_coordsets(RMSDCalculatorType type, double precision_check){
	print_test_tittle(__FUNCTION__);
	cout<<"- Using "<<calculatorTypeToString(type)<<" (prec. "<<precision_check<<"):"<<endl;

	int number_of_coordsets = 5;
	int number_of_atoms = 3239;
	int number_of_CAs = 224;
	// rmsds 1 vs others [ 0.        ,  1.86745533,  2.07960877,  3.60601759,  2.18942902]
	vector<double> not_aligned_CA;
	vector<double> not_aligned_coordinates;
	vector<double> aligned_coordinates;
	double rmsds[5];
	RMSDTools::initializeTo(rmsds, 0., 5);

	load_vector(not_aligned_CA, "data/ligand_mini_CAs");
	load_vector(not_aligned_coordinates, "data/ligand_mini_all");
	load_vector(aligned_coordinates, "data/ligand_mini_all_aligned");

	test_vector_len(not_aligned_CA, number_of_CAs*number_of_coordsets*3, "not aligned CAs");
	test_vector_len(not_aligned_coordinates, number_of_atoms*number_of_coordsets*3, "all not aligned atoms");
	test_vector_len(aligned_coordinates, number_of_atoms*number_of_coordsets*3, "all aligned atoms");

	// RMSD of all atoms using CA for superposition
	RMSDCalculator* calculator1 = RMSDCalculatorFactory::createCalculator(
		type,
		number_of_coordsets,
		number_of_CAs,
		&(not_aligned_CA[0]));
	calculator1->setCalculationCoordinates(number_of_atoms, &(not_aligned_coordinates[0]));
	calculator1->oneVsFollowing(0, rmsds);
	//print_vector("rmsd:",rmsds, 5);
	double expected_rmsds []= {1.864003731005552, 2.076760850428891, 3.596135117728627, 2.182685209336899, 0};
	compareVectors("\tTesting RMSD 1: ", expected_rmsds, rmsds, number_of_coordsets, precision_check); // Only the fourth decimal, as it was obtained with cout without more precission:P
	delete calculator1;

	// RMSD of CA using CA for superposition (default behavior)
	RMSDCalculator* calculator2 = RMSDCalculatorFactory::createCalculator(
		type,
		number_of_coordsets,
		number_of_CAs,
		&(not_aligned_CA[0]));
	calculator2->oneVsFollowing(0, rmsds);
//	print_vector("rmsd:",rmsds, 5);
	double expected_rmsds_2 []= {0.767947519172927, 0.8838644164683896, 0.4177715823462121, 0.3383320758562839, 0};
	compareVectors("\tTesting RMSD 2: ", expected_rmsds_2, rmsds, number_of_coordsets, precision_check);
	delete calculator2;

	// RMSD  of CA using CA for superposition (using the same selection and RMSD subsets)
	RMSDCalculator* calculator3 = RMSDCalculatorFactory::createCalculator(
		type,
		number_of_coordsets,
		number_of_CAs,
		&(not_aligned_CA[0]));
	calculator3->setCalculationCoordinates(number_of_CAs, &(not_aligned_CA[0]));
	calculator3->oneVsFollowing(0, rmsds);
//	print_vector("rmsd:",rmsds, 5);
	compareVectors("\tTesting RMSD 3: ", expected_rmsds_2, rmsds, number_of_coordsets, precision_check);
	delete calculator3;
}

// Fine grain test of qcp with data from the original files in http://theobald.brandeis.edu/qcp/
void test_QCP_Kernel(){
	print_test_tittle(__FUNCTION__);

	int atoms_len = 7;

	double frag_a [] =  {-2.803, -15.373, 24.556, 0.893, -16.062, 25.147,  1.368, -12.371, 25.885, -1.651, -12.153, 28.177, -0.440,
	-15.218, 30.068,  2.551, -13.273, 31.372,  0.105, -11.330, 33.567};

	double frag_b [] =  {-14.739, -18.673,  15.040, -12.473, -15.810,  16.074, -14.802, -13.307,  14.408, -17.782, -14.852,  16.171, -16.124, -14.617,
	19.584, -15.029, -11.037,  18.902, -18.577, -10.001,  17.996};

	double frag_b_copy[atoms_len*3];

	double expected_rotation [] = {0.7221635837820651, 0.6911893731786908, -0.02714790348982324,
								-0.5203825657891069, 0.5170083254696894, -0.6796354733368274,
								-0.4557211246823982, 0.5049352847641727, 0.7330474846272469};

	double expected_rotated_coordinates [] = {
			-2.495176411277905, -1.614342696222044, -4.102116562817358,
			1.092050512774367, -2.016077833910722, -2.931179811963272,
			1.185406934426247, 1.622237699041897, -1.827209404202237,
			-2.082389880657943, 1.176002542749938, 0.04307724778849886,
			-0.8152689506610532, -1.884890665341617, 1.908042480017456,
			2.46847299974008, -0.1403083768834845, 2.716757783430188,
			0.6469047956562158, 2.857379330566034, 4.192628267746734
	};

	QCPSerialKernel kernel;
	double rot_matrix [9];
	double translations[3];


	// Do it step by step
	RMSDTools::copyArrays(frag_b_copy,frag_b,atoms_len*3);
	RMSDTools::centerAllAtOrigin(atoms_len,1,frag_a,translations);
	RMSDTools::centerAllAtOrigin(atoms_len,1,frag_b_copy,translations);

	double rmsd = kernel.calcRMSDOfTwoConformations(frag_a,frag_b_copy,atoms_len,rot_matrix);
	double expected_rmsd =  0.719106;
	compareVectors("\tTesting RMSD: ", &expected_rmsd, &rmsd, 1, 1e-6);

	compareVectors("\tTesting rotation matrix: ", expected_rotation, rot_matrix, 9, 1e-14);

	RMSDTools::rotate3D(atoms_len, frag_b_copy, rot_matrix);
	compareVectors("\tTesting rotated coordinates: ", expected_rotated_coordinates, frag_b_copy, atoms_len*3, 1e-14);

	// Using the function modifying coords
	RMSDTools::copyArrays(frag_b_copy,frag_b,atoms_len*3);
	RMSDTools::centerAllAtOrigin(atoms_len,1,frag_b_copy,translations);
	kernel.oneVsFollowingFitEqualCalcWithConfRotation(
			frag_a,
			-1,
			&rmsd,
			1,
			atoms_len*3,
			atoms_len,
			frag_b_copy);
	compareVectors("\tTesting rotated coordinates: ", expected_rotated_coordinates, frag_b_copy, atoms_len*3, 1e-14);
}

void test_KABSCH_Kernel(){
	print_test_tittle(__FUNCTION__);

		int atoms_len = 7;

		double frag_a [] =  {-2.803, -15.373, 24.556, 0.893, -16.062, 25.147,  1.368, -12.371, 25.885, -1.651, -12.153, 28.177, -0.440,
		-15.218, 30.068,  2.551, -13.273, 31.372,  0.105, -11.330, 33.567};

		double frag_b [] =  {-14.739, -18.673,  15.040, -12.473, -15.810,  16.074, -14.802, -13.307,  14.408, -17.782, -14.852,  16.171, -16.124, -14.617,
		19.584, -15.029, -11.037,  18.902, -18.577, -10.001,  17.996};

		double frag_b_copy[atoms_len*3];

		double expected_rotation [] = {0.7221635837820651, 0.6911893731786908, -0.02714790348982324,
									-0.5203825657891069, 0.5170083254696894, -0.6796354733368274,
									-0.4557211246823982, 0.5049352847641727, 0.7330474846272469};

		double expected_rotated_coordinates [] = {
				-2.495176411277905, -1.614342696222044, -4.102116562817358,
				1.092050512774367, -2.016077833910722, -2.931179811963272,
				1.185406934426247, 1.622237699041897, -1.827209404202237,
				-2.082389880657943, 1.176002542749938, 0.04307724778849886,
				-0.8152689506610532, -1.884890665341617, 1.908042480017456,
				2.46847299974008, -0.1403083768834845, 2.716757783430188,
				0.6469047956562158, 2.857379330566034, 4.192628267746734
		};

		KABSCHSerialKernel kernel;
		double rot_matrix [9];
		double translations[3];
		double U[3][3];

		// Do it step by step
		RMSDTools::copyArrays(frag_b_copy,frag_b,atoms_len*3);
		RMSDTools::centerAllAtOrigin(atoms_len,1,frag_a,translations);
		RMSDTools::centerAllAtOrigin(atoms_len,1,frag_b_copy,translations);

		double rmsd = kernel.calculate_rotation_rmsd(frag_a,frag_b_copy,atoms_len,U);

		double expected_rmsd =  0.719106;
		compareVectors("\tTesting RMSD: ", &expected_rmsd, &rmsd, 1, 1e-6);

		compareVectors("\tTesting rotation matrix: ", expected_rotation, rot_matrix, 9, 1e-14);

		RMSDTools::rotate3D(atoms_len, frag_b_copy, rot_matrix);
		compareVectors("\tTesting rotated coordinates: ", expected_rotated_coordinates, frag_b_copy, atoms_len*3, 1e-14);

		// Using the function modifying coords
		RMSDTools::copyArrays(frag_b_copy,frag_b,atoms_len*3);
		RMSDTools::centerAllAtOrigin(atoms_len,1,frag_b_copy,translations);
		kernel.oneVsFollowingFitEqualCalcWithConfRotation(
				frag_a,
				-1,
				&rmsd,
				1,
				atoms_len*3,
				atoms_len,
				frag_b_copy);
		compareVectors("\tTesting rotated coordinates: ", expected_rotated_coordinates, frag_b_copy, atoms_len*3, 1e-14);

}
