#include "NOSUPOmpKernel.h"
#include "../RMSDCalculationData.h"
#include "../RMSDTools.h"
#include <omp.h>

#include <iostream>
using namespace std;

NoSuperpositionOmpKernel::NoSuperpositionOmpKernel(int number_of_threads){
	this->number_of_threads = number_of_threads;
	omp_set_num_threads(number_of_threads);
}

NoSuperpositionOmpKernel::~NoSuperpositionOmpKernel(){}

void NoSuperpositionOmpKernel::oneVsFollowingFitEqualCalcCoords(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double* second_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(reference,
					second_conformation_coords,
					data->atomsPerFittingConformation);
		}
	}
}

void NoSuperpositionOmpKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	#pragma omp parallel for shared(rmsd)
	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double* calc_conformation_coords =  data->getCalculationConformationAt(second_conformation_index);
		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}
