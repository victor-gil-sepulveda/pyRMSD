/*
 * QTRFITOmpKernel.h
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#ifndef QTRFITOMPKERNEL_H_
#define QTRFITOMPKERNEL_H_

#include "QTRFITSerialKernel.h"

class QTRFITOmpKernel: public QTRFITSerialKernel {
	public:
		QTRFITOmpKernel();
		virtual ~QTRFITOmpKernel();

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

#endif /* QTRFITOMPKERNEL_H_ */
