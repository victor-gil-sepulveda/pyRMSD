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

		virtual void changeCalculationCoords(
				double* calcCoords,
				int number_of_atoms,
				int numberOfConformations){}

		virtual void oneVsFollowingFitEqualCalcWithoutConfRotation(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates) = 0;

		virtual void oneVsFollowingFitEqualCalcWithConfRotation(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates) = 0;

		virtual void oneVsFollowingFitDiffersCalcWithoutConfRotation(
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

		virtual void oneVsFollowingFitDiffersCalcWithConfRotation(
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

		virtual inline void matrixEnd(double* rmsds_tmp, int rmsds_tmp_len, std::vector<double>& rmsds){
			for (int i = 0; i < rmsds_tmp_len; ++i){
				rmsds.push_back(rmsds_tmp[i]);
			}
		}

		virtual void matrixOneVsFollowingFitEqualCalcWithoutConfRotation(
													double* reference,
													int reference_conformation_number,
													double* rmsd,
													int numberOfConformations,
													int coordinatesPerConformation,
													int atomsPerConformation,
													double *allCoordinates){

			oneVsFollowingFitEqualCalcWithoutConfRotation(
									reference,
									reference_conformation_number,
									rmsd,
									numberOfConformations,
									coordinatesPerConformation,
									atomsPerConformation,
									allCoordinates);
		}

		virtual void matrixOneVsFollowingFitDiffersCalcWithoutConfRotation(
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
			oneVsFollowingFitDiffersCalcWithoutConfRotation(
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
