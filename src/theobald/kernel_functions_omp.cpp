#include "kernel_functions_serial.h"
#include "kernel_functions_omp.h"

#include <cmath>
#include <iostream>

using namespace std;

#define floating_point_type double

void ThRMSDSerialOmpKernel::calcRMSDOfOneVsOthers(floating_point_type* all_coordinates,
									 int base_conformation_id,
									 int other_conformations_starting_id,
									 int number_of_conformations,
									 int number_of_atoms,
									 floating_point_type* rmsd){

	int coordinates_per_conformation = number_of_atoms * 3;
	floating_point_type* first_conformation_coords = &(all_coordinates[base_conformation_id*coordinates_per_conformation]);

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_id = other_conformations_starting_id;
			second_conformation_id < number_of_conformations; ++second_conformation_id){
		floating_point_type* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);
		rmsd[second_conformation_id] = ThRMSDSerialKernel::calcRMSDOfTwoConformations(first_conformation_coords, second_conformation_coords, number_of_atoms);
	}
}
