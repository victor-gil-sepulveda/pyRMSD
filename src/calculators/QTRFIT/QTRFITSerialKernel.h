/*
 * QTRFITSerialKernel.h
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#ifndef QTRFITSERIALKERNEL_H_
#define QTRFITSERIALKERNEL_H_

#include "../KernelFunctions.h"

class QTRFITSerialKernel: public KernelFunctions {
	public:
		QTRFITSerialKernel();
		virtual ~QTRFITSerialKernel();

		virtual void oneVsFollowingFitEqualCalcWithoutConfRotation(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates);

		virtual void oneVsFollowingFitEqualCalcWithConfRotation(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				int numberOfConformations,
				int coordinatesPerConformation,
				int atomsPerConformation,
				double *allCoordinates);

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
				double *allRMSDCoordinates);

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
				double *allRMSDCoordinates);

};

#endif /* QTRFITSERIALKERNEL_H_ */
