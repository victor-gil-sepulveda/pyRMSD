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

#define TODOUBLEP(vec) (&(vec[0]))

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

void test_mean_coordinates(){
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



void test_center_coordinates(){
	print_test_tittle(__FUNCTION__);

	vector<double> coordinates, centered_coordinates;
	vector<int> shape, centered_shape;
	double calculated_centers[18];
	double expected_centers [] = {
			  14.3707713,   47.34880717,  25.46220179,
			   7.3707713,    2.3488296,   31.46225112,
			 -10.62925112,  -6.65117489,  30.46225112,
			  15.37077578,  21.3488565,   46.46220179,
			  -5.62919731,   3.34880269,  21.46222422,
			   8.37079821,  96.34881166,  25.46223767};

	RMSDTools::initializeTo(calculated_centers, 0, 18);

	// Not centered coordinates
	load_pdb_coords(	coordinates,
						shape,
						"test_data/Models/prot_plus_ligand_similar/prot_plus_ligand_offset.CA.coords");

	for (int i = 0; i < shape[0]; ++i){
		RMSDTools::geometricCenter(shape[1], TODOUBLEP(coordinates)+(i*shape[1]*3) , calculated_centers+(i*3));
	}

	// Centers must be equal to the expected ones
	compareVectors("\tGeometric centers are as expected: ",
					expected_centers,
					calculated_centers,
					18,
					1e-7);


	// Load the coordinates centered with Python and center them
	load_pdb_coords(centered_coordinates,
						centered_shape,
						"test_data/Models/prot_plus_ligand_similar/prot_plus_ligand_offset.CA.centered.coords");

	RMSDTools::centerAllAtOrigin(shape[1],shape[0],TODOUBLEP(coordinates));

	// Coordinates must coincide
	compareVectors("\tCoordinates have been centered: ",
			TODOUBLEP(centered_coordinates),
			TODOUBLEP(coordinates),
			shape[0]*shape[1]*shape[2],
			1e-12);


}

void test_superposition_with_fit(	RMSDCalculatorType type,
										const char* initial_coords_file,
										const char* final_coords_file,
										const char* rmsd_results_file,
										double precision_of_check){
	print_test_tittle(__FUNCTION__);
	print_calculator_and_precission(type, precision_of_check);

	vector<double> not_superposed_fit_coordinates,
					expected_superposed_fit_coordinates,
					calculated_rmsds, expected_rmsds;

	vector<int> expected_superposed_fit_coordinates_shape,
				not_superposed_fit_coordinates_shape;

	load_vector(expected_rmsds, rmsd_results_file);

	load_pdb_coords(not_superposed_fit_coordinates,
						not_superposed_fit_coordinates_shape,
						initial_coords_file);

	// Prody's results are superposed but the centering has been canceled,
	// it is necessary then to move then again to their original places
	load_and_center_pdb_coords(expected_superposed_fit_coordinates,
						expected_superposed_fit_coordinates_shape,
						final_coords_file);

	calculated_rmsds.resize(not_superposed_fit_coordinates_shape[0],0);
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
									type,
									not_superposed_fit_coordinates_shape[0],
									not_superposed_fit_coordinates_shape[1],
									TODOUBLEP(not_superposed_fit_coordinates));
	calculator->oneVsFollowing(0, TODOUBLEP(calculated_rmsds));

	// RMSDs must be the same
	compareVectors("\tCalculated RMSDs coincide with golden: ",
			&(expected_rmsds[1]),
			TODOUBLEP(calculated_rmsds),
			not_superposed_fit_coordinates_shape[0]-1, precision_of_check);

	// Final centered coordinates must be the superposed coordinates
	compareVectors("\tInitial coordinates have been superposed: ",
			TODOUBLEP(expected_superposed_fit_coordinates),
			TODOUBLEP(not_superposed_fit_coordinates),
				not_superposed_fit_coordinates_shape[0] *
				not_superposed_fit_coordinates_shape[1] *
				not_superposed_fit_coordinates_shape[2],
			precision_of_check);

	//save_pdb_coords(not_superposed_fit_coordinates,not_superposed_fit_coordinates_shape,"calculated.coords");
	//save_pdb_coords(expected_superposed_fit_coordinates,expected_superposed_fit_coordinates_shape,"expected.coords");
	delete calculator;
}


void test_superposition_with_fit_and_calc(RMSDCalculatorType type,
												const char* initial_prot_coords_file,
												const char* final_prot_coords_file,
												const char* initial_lig_coords_file,
												const char* final_lig_coords_file,
												const char* rmsd_results_file,
												double precision_of_check){
	print_test_tittle(__FUNCTION__);
	print_calculator_and_precission(type, precision_of_check);

	vector<double>     not_superposed_fit_coordinates,
						expected_superposed_fit_coordinates,
						not_superposed_calc_coordinates,
						expected_superposed_calc_coordinates,
						calculated_rmsds, expected_rmsds,
						centers;

	vector<int> expected_superposed_fit_coordinates_shape,
				not_superposed_fit_coordinates_shape,
				expected_superposed_calc_coordinates_shape,
				not_superposed_calc_coordinates_shape;

	load_vector(expected_rmsds, rmsd_results_file);

	load_pdb_coords(not_superposed_fit_coordinates,
						not_superposed_fit_coordinates_shape,
						initial_prot_coords_file);

	// Prody's results are superposed but the centering has been canceled,
	// it is necessary then to move then again to their original places
	load_and_center_pdb_coords(expected_superposed_fit_coordinates,
						expected_superposed_fit_coordinates_shape,
						final_prot_coords_file,
						&centers);

	load_pdb_coords(not_superposed_calc_coordinates,
						not_superposed_calc_coordinates_shape,
						initial_lig_coords_file);

	// The case of a different non-centered calculation coordset is a little bit more tricky,
	// to preserve the relative distance to the center, one has to move this coordinates using
	// the same centers got for the fitting coordinates
	load_and_move_pdb_coords(expected_superposed_calc_coordinates,
						expected_superposed_calc_coordinates_shape,
						final_lig_coords_file,
						TODOUBLEP(centers));

	calculated_rmsds.resize(not_superposed_fit_coordinates_shape[0],0);
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
									type,
									not_superposed_fit_coordinates_shape[0],
									not_superposed_fit_coordinates_shape[1],
									TODOUBLEP(not_superposed_fit_coordinates));
	calculator->setCalculationCoordinates(
			not_superposed_calc_coordinates_shape[1],
			TODOUBLEP(not_superposed_calc_coordinates));

	calculator->oneVsFollowing(0, TODOUBLEP(calculated_rmsds));

//	print_vector("expected RMSD: ", TODOUBLEP(expected_rmsds),expected_rmsds.size(),8);
//	print_vector("calcted. RMSD: ", TODOUBLEP(calculated_rmsds), calculated_rmsds.size(),8);

	// RMSDs must be the same
	compareVectors("\tCalculated RMSDs coincide with golden: ",
			&(expected_rmsds[1]),
			TODOUBLEP(calculated_rmsds),
			not_superposed_fit_coordinates_shape[0]-1, precision_of_check);

	compareVectors("\tInitial fitting coordinates have been superposed: ",
			TODOUBLEP(expected_superposed_fit_coordinates),
			TODOUBLEP(not_superposed_fit_coordinates),
				not_superposed_fit_coordinates_shape[0] *
				not_superposed_fit_coordinates_shape[1] *
				not_superposed_fit_coordinates_shape[2],
			precision_of_check);

	compareVectors("\tAnd also calculation coordinates: ",
			TODOUBLEP(expected_superposed_calc_coordinates),
			TODOUBLEP(not_superposed_calc_coordinates),
				not_superposed_calc_coordinates_shape[0] *
				not_superposed_calc_coordinates_shape[1] *
				not_superposed_calc_coordinates_shape[2],
			precision_of_check);
}

void test_iterative_superposition_with_fit(RMSDCalculatorType type,
				const char* initial_prot_coords_file,
				const char* final_prot_coords_file,
				const char* iteration_rmsd_results_file,
				double precision_of_check){

	print_test_tittle(__FUNCTION__);
	print_calculator_and_precission(type, precision_of_check);

	vector<double> 		not_iteratively_superposed_fit_coordinates,
						expected_iteratively_superposed_fit_coordinates,
						calculated_iteration_rmsds, expected_iteration_rmsds,
						iterposition_step;

	double* reference_coords = NULL;
	double* mean_coords = NULL;

	vector<int> expected_fit_coordinates_shape,
				fit_coordinates_shape,
				one_step_fit_size;

	load_vector(expected_iteration_rmsds, iteration_rmsd_results_file);

	load_pdb_coords(not_iteratively_superposed_fit_coordinates,
						fit_coordinates_shape,
						initial_prot_coords_file);

	// Prody's results are superposed but the centering has been canceled,
	// it is necessary then to move then again to their original places
	load_and_center_pdb_coords(expected_iteratively_superposed_fit_coordinates,
								expected_fit_coordinates_shape,
								final_prot_coords_file);

	// Step by step test
	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
									type,
									fit_coordinates_shape[0],
									fit_coordinates_shape[1],
									TODOUBLEP(not_iteratively_superposed_fit_coordinates));

	reference_coords = new double[fit_coordinates_shape[1]*3];
	RMSDTools::initializeTo(reference_coords,0.,fit_coordinates_shape[1]*3);
	mean_coords = new double[fit_coordinates_shape[1]*3];
	RMSDTools::initializeTo(mean_coords,0.,fit_coordinates_shape[1]*3);

	RMSDTools::copyArrays(reference_coords,
			TODOUBLEP(not_iteratively_superposed_fit_coordinates),
			fit_coordinates_shape[1]*3);

	calculator->iterativeSuperpositionStep(reference_coords, mean_coords);
	//sudo mousepad /usr/local/lib/python2.7/dist-packages/prody/ensemble/ensemble.py
	load_and_center_pdb_coords(iterposition_step,
							   one_step_fit_size,
							   "test_data/Iterpose_Fit_CA_Rot_BEN/steps/stretching_trajectory_offset_ligand.iterposed_CA.1_step.coords");

	compareVectors("\tFitting coordinates and golden are equal: ",
					TODOUBLEP(iterposition_step),
					TODOUBLEP(not_iteratively_superposed_fit_coordinates),
						fit_coordinates_shape[0] *
						fit_coordinates_shape[1] *
						fit_coordinates_shape[2],
					precision_of_check);

	compareVectors("\tMean and reference is the same for step 0: ",
						mean_coords,
						reference_coords,
							one_step_fit_size[1] *
							one_step_fit_size[2],
						precision_of_check);

	load_and_center_pdb_coords(iterposition_step,
								   one_step_fit_size,
								   "test_data/Iterpose_Fit_CA_Rot_BEN/mean/mean_step0.coords");
	compareVectors("\tMean is the correct mean and reference is the new mean for step 0: ",
					TODOUBLEP(iterposition_step),
					reference_coords,
						one_step_fit_size[1] *
						one_step_fit_size[2],
					precision_of_check);
	delete [] reference_coords;
	delete [] mean_coords;

//	calculated_iteration_rmsds.resize(200,0);
//	RMSDCalculator* calculator = RMSDCalculatorFactory::createCalculator(
//									type,
//									fit_coordinates_shape[0],
//									fit_coordinates_shape[1],
//									TODOUBLEP(not_iteratively_superposed_fit_coordinates));
//
//	calculator->iterativeSuperposition(1e-4, TODOUBLEP(calculated_iteration_rmsds));
//
//	print_vector("calculated RMSD: ", TODOUBLEP(calculated_iteration_rmsds), calculated_iteration_rmsds.size(),8);
//	print_vector("expected RMSD: ", TODOUBLEP(expected_iteration_rmsds), expected_iteration_rmsds.size(),8);
//
//	compareVectors("\tFitting coordinates and golden are equal: ",
//				TODOUBLEP(expected_iteratively_superposed_fit_coordinates),
//				TODOUBLEP(not_iteratively_superposed_fit_coordinates),
//					fit_coordinates_shape[0] *
//					fit_coordinates_shape[1] *
//					fit_coordinates_shape[2],
//				precision_of_check);

}
