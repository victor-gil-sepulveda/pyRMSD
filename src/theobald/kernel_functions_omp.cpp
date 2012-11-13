#include "kernel_functions_serial.h"
#include "kernel_functions_omp.h"

#include <cmath>
#include <iostream>

using namespace std;

#define floating_point_type double

///////////////////////////////////////////////////////////////
/// \remarks
/// OpenMP naive port of the serial version for the same function. Parallelism is at 'calcRMSDOfTwoConformations'. This
/// approach gives quite decent speedups anyway (and surely a very good work load).
/// Given a trajectory with N conformations stored in 'all_coordinates', it calculates the rmsd between the conformation
/// with 'i' = 'base_conformation_id' and the conformations in the range [j+1,M)
///
///
/// \param 	all_coordinates [In] The array storing all coordinates for all conformations in order:
/// (conf0_atom_0_x,conf0_atom_0_y,conf0_atom_0_z,conf0_atom_1_x,conf0_atom_1_y ... confN-1_atom_N_z)
///
/// \param 	base_conformation_id [In] The reference conformation id ('i' in the above explanation).
///
/// \param 	other_conformations_starting_id [In] The value of j in the above explanation.
///
/// \param 	number_of_conformations [In] Number of conformations stored in the coordinates array.
///
/// \param 	number_of_atoms [In] Number of atoms of any of the conformations.
///
/// \param 	rmsd [In] Array where the rmsds are stored (rmsd[j]  is the rmsd between conformation 'i' and
///	conformation 'j').
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
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
