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


		virtual void superpose(unsigned int n, double * const coord_fit, double* const coord_ref);

		virtual void superpose(unsigned int fit_n, double * const fit_coords, double* const fit_ref_coords,
				unsigned int rmsd_n, double* const calc_coords);

		virtual void superpositionQuatFit(unsigned int n, const double * const x, const double * const y, double q[4], double u[3][3]);

		virtual void generateUpperQuadraticMatrix(double upperQuadMatrix[4][4], unsigned int n, const double * const y, const double *const x);

		virtual void computeFittedQuaternionFromUpperQuadraticMatrix(double upperQuadMatrix[4][4], double q[4]);

		virtual void generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3]);

};

#endif /* QTRFITSERIALKERNEL_H_ */
