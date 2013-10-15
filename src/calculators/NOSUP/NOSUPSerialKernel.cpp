#include "NOSUPSerialKernel.h"
#include "../RMSDCalculationData.h"
#include "../RMSDTools.h"

/*
 * Implementation for cases in which superposition is not required (structures are already superposed).
 * In this cases superposition may converge very fast, but it will still add an important overhead. Indeed
 * the calculator class will center structures when necessary, adding an overhead (which will be minimum when
 * calculating the RMSD matrix).
 * Due to the centering step 'oneVsFollowingFitDiffersCalcCoords' needs to have the two coordinate sets defined
 * or the result will be undetermined.
 */

NoSuperpositionSerialKernel::NoSuperpositionSerialKernel(){}

NoSuperpositionSerialKernel::~NoSuperpositionSerialKernel(){}

void NoSuperpositionSerialKernel::oneVsFollowingFitEqualCalcCoords(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			RMSDCalculationData* data){

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

void NoSuperpositionSerialKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

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
