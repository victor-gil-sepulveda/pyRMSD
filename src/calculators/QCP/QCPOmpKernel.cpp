#include <cmath>
#include <iostream>
#include "QCPOmpKernel.h"
#include "../RMSDTools.h"
#include <omp.h>
#include "../RMSDCalculationData.h"
using namespace std;

QCPOmpKernel::QCPOmpKernel(int number_of_threads){
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Given a trajectory with N conformations stored in 'all_coordinates', it calculates the rmsd between the conformation
/// with 'i' = 'base_conformation_id' and the conformations in the range [j+1,M)
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
void QCPOmpKernel::oneVsFollowingFitEqualCalcCoords(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){
		double rot_matrix[9];

		double* second_conformation_coords = data->getFittingConformationAt(second_conformation_index);

		RMSDTools::initializeTo(rot_matrix, 0.0, 9);

		double rmsd_val = calcRMSDOfTwoConformations(	reference,
													  	second_conformation_coords,
														data->atomsPerFittingConformation,
														rot_matrix);
		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = rmsd_val;
		}

		RMSDTools::rotate3D(data->atomsPerFittingConformation,
				second_conformation_coords,
				rot_matrix);
	}
}


void QCPOmpKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double rot_matrix[9];

		double* fit_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		double* calc_conformation_coords =  data->getCalculationConformationAt(second_conformation_index);

		RMSDTools::initializeTo(rot_matrix, 0.0, 9);

		calcRMSDOfTwoConformations(	fitReference,
									fit_conformation_coords,
									data->atomsPerFittingConformation,
									rot_matrix);


		RMSDTools::rotate3D(data->atomsPerFittingConformation,
				fit_conformation_coords,
				rot_matrix);

		RMSDTools::rotate3D(data->atomsPerCalculationConformation,
				calc_conformation_coords,
				rot_matrix);

		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}
