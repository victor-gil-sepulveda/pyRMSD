/*
 * KernelFunctions.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef KERNELFUNCTIONS_H_
#define KERNELFUNCTIONS_H_

#include <vector>

class KernelFunctions{
	public:
		KernelFunctions(){}
		virtual ~KernelFunctions(){}

		virtual void setCalculationCoords(
				double* calcCoords,
				int number_of_atoms,
				int numberOfConformations){}

		virtual void oneVsFollowingFitEqualCalcCoords(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates) = 0;

		virtual void oneVsFollowingFitDiffersCalcCoords(
				double* fitReference,
				double* calcReference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates,
				int coordinatesPerRMSDConformation,
				int atomsPerRMSDConformation,
				double *allRMSDCoordinates) = 0;

		virtual void matrixInit(
				double* allFittingCoordinates,
				int coordinatesPerFittingConformation,
				double* allCalculationCoordinates,
				int coordinatesPerCalculationConformation,
				int numberOfConformations){

		}

		virtual void matrixEnd(
						int rmsds_tmp_len,
						std::vector<double>& rmsds){
		}

		virtual void matrixOneVsFollowingFitEqualCalc(
													double* reference,
													int reference_conformation_number,
													double* rmsd,
													int numberOfConformations,
													int coordinatesPerConformation,
													int atomsPerConformation,
													double *allCoordinates){

			oneVsFollowingFitEqualCalcCoords(
									reference,
									reference_conformation_number,
									rmsd,
									numberOfConformations,
									coordinatesPerConformation,
									atomsPerConformation,
									allCoordinates);
		}

		virtual void matrixOneVsFollowingFitDiffersCalc(
														double* fitReference,
														double* calcReference,
														int reference_conformation_number,
														double* rmsd,
														int numberOfConformations,
														int coordinatesPerConformation,
														int atomsPerConformation,
														double *allCoordinates,
														int coordinatesPerRMSDConformation,
														int atomsPerRMSDConformation,
														double *allRMSDCoordinates){
			oneVsFollowingFitDiffersCalcCoords(
												fitReference,
												calcReference,
												reference_conformation_number,
												rmsd,
												numberOfConformations,
												coordinatesPerConformation,
												atomsPerConformation,
												allCoordinates,
												coordinatesPerRMSDConformation,
												atomsPerRMSDConformation,
												allRMSDCoordinates);
		}
};


#endif /* KERNELFUNCTIONS_H_ */
