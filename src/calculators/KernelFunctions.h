/*
 * KernelFunctions.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef KERNELFUNCTIONS_H_
#define KERNELFUNCTIONS_H_

class KernelFunctions{
	public:
		KernelFunctions(){}
		virtual ~KernelFunctions(){}

		virtual void init(
						double* coordinates,
						int atomsPerConformation,
						int coordinatesPerConformation,
						int numberOfConformations){}

		virtual void changeCalculationCoords(double*){}

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

		virtual void oneVsAllFitDiffersCalcWithConfRotation(
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
};


#endif /* KERNELFUNCTIONS_H_ */
