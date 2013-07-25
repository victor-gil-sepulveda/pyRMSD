/*
 * QTRFITSerialKernel.cpp
 *
 *  Created on: 05/03/2013
 *      Author: victor
 */

#include "QTRFITSerialKernel.h"
#include "../RMSDTools.h"
#include <cstddef>
#include <iostream>
#include "../RMSDCalculationData.h"

using namespace std;

QTRFITSerialKernel::QTRFITSerialKernel() {}

QTRFITSerialKernel::~QTRFITSerialKernel() {}

void QTRFITSerialKernel::oneVsFollowingFitEqualCalcCoords(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data) {

	for (int i = reference_conformation_number + 1; i < data->numberOfConformations; ++i) {
		// Real conformation coordinates are used, so they are modified
		double* conformation_coords = data->getFittingConformationAt(i);

		superpose(data->atomsPerFittingConformation,
				conformation_coords,
				reference);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL){
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					reference,
					conformation_coords,
					data->atomsPerFittingConformation);
		}
	}
}

void QTRFITSerialKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data) {

	for (int i = reference_conformation_number + 1; i < data->numberOfConformations;	++i) {

		double* fit_conformation_coords = data->getFittingConformationAt(i);
		double* calc_conformation_coords = data->getCalculationConformationAt(i);

		superpose(data->atomsPerFittingConformation,
				fit_conformation_coords,
				fitReference,
				data->atomsPerCalculationConformation,
				calc_conformation_coords);

		// rmsd vector can be null if we are only interested in conformation superposition
		if (rmsd != NULL) {
			rmsd[i - (reference_conformation_number + 1)] = RMSDTools::calcRMS(
					calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}

// After the function both coordinate sets are modified
void QTRFITSerialKernel::superpose(unsigned int n, double * const coord_fit, double* const coord_ref){
	double rot_matrix[3][3];
	double quaternion[4];

	// Get rotation matrix
	superpositionQuatFit(n, coord_fit, coord_ref, quaternion, rot_matrix);

	// Rotate the reference molecule by the rotation matrix u
	RMSDTools::rotate3D(n, coord_fit, rot_matrix);

}

// superposes calc_coords over ref_coords and leaves them @ origin
// Once finished,
void QTRFITSerialKernel::superpose(unsigned int fit_n, double * const fit_coords, double* const fit_ref_coords,
							unsigned int calc_n, double* const calc_coords){
	double rot_matrix[3][3];
	double quaternion[4];

	// Get rotation matrix
	superpositionQuatFit(fit_n,  fit_coords, fit_ref_coords, quaternion, rot_matrix);

	// Rotate the calculation and fit coordinate sets by using rotation matrix u
	RMSDTools::rotate3D(fit_n, fit_coords, rot_matrix);
	RMSDTools::rotate3D(calc_n, calc_coords, rot_matrix);
}

void QTRFITSerialKernel::superpositionQuatFit(unsigned int n, const double * const structure_to_fit, const double * const reference, double q[4], double u[3][3])
{
	//vector<double> upperQuadMatrix(16, 0);
	double upperQuadMatrix[4][4];

	// Generate the upper triangle of the quadratic matrix
	generateUpperQuadraticMatrix(upperQuadMatrix, n, reference, structure_to_fit);

	// Compute quaternion
	computeFittedQuaternionFromUpperQuadraticMatrix(&upperQuadMatrix[0], q);

	// Generate the rotation matrix
	generateLeftRotationMatrixFromNormalizedQuaternion(q, u);
}

void QTRFITSerialKernel::generateUpperQuadraticMatrix(double upperQuadMatrix[4][4], unsigned int n, const double * const reference, const double *const structure_to_fit){

	const double* y = reference;
	const double* const x = structure_to_fit;

	//Initialize upperQuadMatrix
	for(unsigned int i = 0; i < 4; ++i){
		for(unsigned int j = 0; j < 4; ++j){
			upperQuadMatrix[i][j] = 0.0;
		}
	}

	// Generate the upper triangle of the quadratic matrix
	double xxyx = 0.0;
	double xxyy = 0.0;
	double xxyz = 0.0;
	double xyyx = 0.0;
	double xyyy = 0.0;
	double xyyz = 0.0;
	double xzyx = 0.0;
	double xzyy = 0.0;
	double xzyz = 0.0;

	for(unsigned int i=0; i<n; ++i)
	{
		unsigned int k = 3 * i;

		xxyx += x[k] * y[k];
		xxyy += x[k] * y[k+1];
		xxyz += x[k] * y[k+2];
		xyyx += x[k+1] * y[k];
		xyyy += x[k+1] * y[k+1];
		xyyz += x[k+1] * y[k+2];
		xzyx += x[k+2] * y[k];
		xzyy += x[k+2] * y[k+1];
		xzyz += x[k+2] * y[k+2];
	}

	upperQuadMatrix[0][0] = xxyx + xyyy + xzyz;

	upperQuadMatrix[0][1] = xzyy - xyyz;
	upperQuadMatrix[1][1] = xxyx - xyyy - xzyz;

	upperQuadMatrix[0][2] = xxyz - xzyx;
	upperQuadMatrix[1][2] = xxyy + xyyx;
	upperQuadMatrix[2][2] = xyyy - xzyz - xxyx;

	upperQuadMatrix[0][3] = xyyx - xxyy;
	upperQuadMatrix[1][3] = xzyx + xxyz;
	upperQuadMatrix[2][3] = xyyz + xzyy;
	upperQuadMatrix[3][3] = xzyz - xxyx - xyyy;
}

void QTRFITSerialKernel::computeFittedQuaternionFromUpperQuadraticMatrix(double upperQuadMatrix[4][4], double q[4]){
	double eigvec[4][4],eigval[4];
	RMSDTools::jacobi(upperQuadMatrix,eigval,eigvec);
	q[0] = eigvec[0][3];
	q[1] = eigvec[1][3];
	q[2] = eigvec[2][3];
	q[3] = eigvec[3][3];
}

void QTRFITSerialKernel::generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3])
{
	u[0][0] = q[0] * q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
	u[0][1] = 2.0 * (q[2] * q[1] + q[0]*q[3]);
	u[0][2] = 2.0 * (q[3] * q[1] - q[0]*q[2]);

	u[1][0] = 2.0 * (q[1] * q[2] - q[0]*q[3]);
	u[1][1] = q[0] * q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
	u[1][2] = 2.0 * (q[3] * q[2] + q[0]*q[1]);

	u[2][0] = 2.0 * (q[1] * q[3] + q[0]*q[2]);
	u[2][1] = 2.0 * (q[2] * q[3] - q[0]*q[1]);
	u[2][2] = q[0] * q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

