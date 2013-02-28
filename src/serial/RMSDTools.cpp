
#include "RMSDTools.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;


// After the function both coordinate sets are modified
void RMSDTools::superpose(unsigned int n, double * const coord_fit, double* const coord_ref){
	//double center_fit[3] = {0.,0.,0.};
	//double center_ref[3] = {0.,0.,0.};
	double u[3][3];
	double q[4];

	// Get rotation matrix
	RMSDTools::superpositionQuatFit(n, coord_fit, coord_ref, q, u);

	// Rotate the reference molecule by the rotation matrix u
	RMSDTools::rotate3D(n, coord_fit, u);

	/*// Shift the reference molecule to the geometric center of the fit
	RMSDTools::shift3D(n, coord_ref, center_ref, +1.);

	// Shift to the origin
	RMSDTools::shift3D(n, coord_fit, center_ref, +1.);*/
}

// superposes calc_coords over ref_coords and leaves them @ origin
// Once finished,
void RMSDTools::superpose(unsigned int fit_n, double * const fit_coords, double* const fit_ref_coords,
							unsigned int calc_n, double* const calc_coords, double* const calc_reference){
	double u[3][3];
	double q[4];

	// Get rotation matrix
	RMSDTools::superpositionQuatFit(fit_n,  fit_coords, fit_ref_coords,q, u);

	// Rotate the calculation coordinate set using rotation matrix u
	RMSDTools::rotate3D(calc_n, calc_coords, u);
}
/*
 * Centers all the conformations to origin and stores the movement vectors.
 *
 */
void RMSDTools::centerAllAtOrigin(unsigned int atomsPerConformation, unsigned int numberOfConformations, double * const all_coords, double* const translations){

	unsigned int coordsPerConformation = atomsPerConformation * 3;

	#pragma omp parallel for
	for (unsigned int i = 0; i < numberOfConformations; ++i){
		double* center = &(translations[3*i]);
		double* coords = &(all_coords[coordsPerConformation*i]);
		RMSDTools::geometricCenter(atomsPerConformation, coords, center);
		RMSDTools::shift3D(atomsPerConformation, coords, center, -1.);
	}
}

void RMSDTools::applyTranslationsToAll(unsigned int atomsPerConformation, unsigned int numberOfConformations, double * const all_coords, double* const translations){
	unsigned int coordsPerConformation = atomsPerConformation * 3;

	#pragma omp parallel for
	for (unsigned int i = 0; i < numberOfConformations; ++i){
		double* translation_vector = &(translations[3*i]);
		double* coords = &(all_coords[coordsPerConformation*i]);
		RMSDTools::shift3D(atomsPerConformation, coords, translation_vector, 1.);
	}
}

void RMSDTools::applyTranslationToAll(unsigned int atomsPerConformation, unsigned int numberOfConformations, double * const all_coords, double* const translation_vector){
	unsigned int coordsPerConformation = atomsPerConformation * 3;

	#pragma omp parallel for
	for (unsigned int i = 0; i < numberOfConformations; ++i){
		double* coords = &(all_coords[coordsPerConformation*i]);
		RMSDTools::shift3D(atomsPerConformation, coords, translation_vector, 1.);
	}
}



void RMSDTools::geometricCenter(unsigned int n, const double * const x, double * const center){
	unsigned int i;

	// Initialize variables before the loop
	center[0] = 0.0;
	center[1] = 0.0;
	center[2] = 0.0;

	for(i=0; i<n; i++)
	{
		int offset = 3*i;
		center[0] += x[offset];
		center[1] += x[offset+1];
		center[2] += x[offset+2];
	}

	center[0] /= n;
	center[1] /= n;
	center[2] /= n;
}

void RMSDTools::shift3D(unsigned int numberOfPoints, double * const x, double trans[3], double scalar){
	double shiftVector[3];

	shiftVector[0] = trans[0]*scalar;
	shiftVector[1] = trans[1]*scalar;
	shiftVector[2] = trans[2]*scalar;

	for(unsigned int i=0; i<numberOfPoints; ++i){
		int offset = 3*i;
		x[offset] += shiftVector[0];
		x[offset+1] += shiftVector[1];
		x[offset+2] += shiftVector[2];
	}
}

void RMSDTools::superpositionQuatFit(unsigned int n, const double * const structure_to_fit, const double * const reference, double q[4], double u[3][3])
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

void RMSDTools::generateUpperQuadraticMatrix(double upperQuadMatrix[4][4], unsigned int n, const double * const reference, const double *const structure_to_fit){

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

void RMSDTools::computeFittedQuaternionFromUpperQuadraticMatrix(double upperQuadMatrix[4][4], double q[4]){
	double eigvec[4][4],eigval[4];
	RMSDTools::jacobi(upperQuadMatrix,eigval,eigvec);
	q[0] = eigvec[0][3];
	q[1] = eigvec[1][3];
	q[2] = eigvec[2][3];
	q[3] = eigvec[3][3];

}

void RMSDTools::generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3])
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


void RMSDTools::rotate3D(unsigned int n, double * const x, double u[3][3]){
	// We go through all selected atoms
	for(unsigned int i=0; i<n; ++i){
		int offset = i*3;
		double x_tmp_0,x_tmp_1,x_tmp_2;
		x_tmp_0 = x[offset];
		x_tmp_1 = x[offset+1];
		x_tmp_2 = x[offset+2];

		// An rotate each of them
		x[offset] 	= u[0][0] * x_tmp_0 + u[0][1] * x_tmp_1 + u[0][2] * x_tmp_2;
		x[offset+1] = u[1][0] * x_tmp_0 + u[1][1] * x_tmp_1 + u[1][2] * x_tmp_2;
		x[offset+2] = u[2][0] * x_tmp_0 + u[2][1] * x_tmp_1 + u[2][2] * x_tmp_2;
	}
}

double RMSDTools::calcRMS(const double * const x, const double * const y, unsigned int num_atoms){
	double sum_res = 0.0;

	for(unsigned int i=0; i<num_atoms*3; ++i){
		sum_res += (x[i] - y[i]) * (x[i] - y[i]);
	}

	return sqrt(sum_res/num_atoms);
}

void RMSDTools::initializeTo(double* array, double value, int array_len){
//	#pragma omp parallel for
//	for(int i = 0; i < array_len; ++i){
//		array[i] = value;
//	}
	fill(array, array+array_len, value);
}

//array1 = array2
void RMSDTools::copyArrays(double* array1, double* array2, int array_len){
//	#pragma omp parallel for
//	for(int i = 0; i < array_len; ++i){
//		array1[i] = array2[i];
//	}
	copy(array2,array2+array_len,array1);
}

void RMSDTools::calculateMeanCoordinates(double* meanCoordinates, double* allCoordinates,
											int numberOfConformations, int atomsPerConformation){

	// Zero mean coordinates
	RMSDTools::initializeTo(meanCoordinates, 0.0, atomsPerConformation*3);

	// Do calculation
	for (int i  = 0; i <  numberOfConformations; ++i){
		int conformation_offset = i*atomsPerConformation*3;
		#pragma omp parallel for shared(meanCoordinates)
		for (int j = 0; j < atomsPerConformation; ++j){
			int atom_offset = 3*j;
			int offset = conformation_offset + atom_offset;
			meanCoordinates[atom_offset] += allCoordinates[ offset ];
			meanCoordinates[atom_offset+1] += allCoordinates[ offset + 1];
			meanCoordinates[atom_offset+2] += allCoordinates[ offset + 2];
		}
	}

	// Divide by the number of conformations
	#pragma omp parallel for shared(meanCoordinates)
	for (int i = 0; i < atomsPerConformation*3; ++i){
		meanCoordinates[i] /= numberOfConformations;
	}
}

/*
 David J. Heisterberg
 The Ohio Supercomputer Center
 1224 Kinnear Rd.
 Columbus, OH  43212-1163
 (614)292-6036
 djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh

 Translated to C from fitest.f program and interfaced with Xmol program
 by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl

 Some minor changes and indentation by Víctor Gil Sepúlveda

 Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
 The program can be copied and distributed freely, provided that
 this copyright in not removed. You may acknowledge the use of the
 program in published material as:
 David J. Heisterberg, 1990, unpublished results.
*/
void RMSDTools::jacobi(double a[4][4], double d[4], double v[4][4], int nrot){
	double onorm, dnorm;
	double b, dma, q, t, c, s;
	double atemp, vtemp, dtemp;
	int i, j, k, l;

	for (j = 0; j <= 3; j++) {
		for (i = 0; i <= 3; i++) {
			v[i][j] = 0.0;
		}
		v[j][j] = 1.0;
		d[j] = a[j][j];
	}

	for (l = 1; l <= nrot; l++) {
		dnorm = 0.0;
		onorm = 0.0;
		for (j = 0; j <= 3; j++) {
			dnorm = dnorm + fabs(d[j]);
			for (i = 0; i <= j - 1; i++) {
				onorm = onorm + fabs(a[i][j]);
			}
		}

		if((onorm/dnorm) <= 1.0e-12){
			break;
		}

		for (j = 1; j <= 3; j++) {
			for (i = 0; i <= j - 1; i++) {
				b = a[i][j];
				if(fabs(b) > 0.0) {
					dma = d[j] - d[i];
					if((fabs(dma) + fabs(b)) <=  fabs(dma)) {
						t = b / dma;
					}
					else {
						q = 0.5 * dma / b;
						t = 1.0/(fabs(q) + sqrt(1.0+q*q));
						if(q < 0.0) {
							t = -t;
						}
					}
					c = 1.0/sqrt(t * t + 1.0);
					s = t * c;
					a[i][j] = 0.0;
					for (k = 0; k <= i-1; k++) {
						atemp = c * a[k][i] - s * a[k][j];
						a[k][j] = s * a[k][i] + c * a[k][j];
						a[k][i] = atemp;
					}
					for (k = i+1; k <= j-1; k++) {
						atemp = c * a[i][k] - s * a[k][j];
						a[k][j] = s * a[i][k] + c * a[k][j];
						a[i][k] = atemp;
					}
					for (k = j+1; k <= 3; k++) {
						atemp = c * a[i][k] - s * a[j][k];
						a[j][k] = s * a[i][k] + c * a[j][k];
						a[i][k] = atemp;
					}
					for (k = 0; k <= 3; k++) {
						vtemp = c * v[k][i] - s * v[k][j];
						v[k][j] = s * v[k][i] + c * v[k][j];
						v[k][i] = vtemp;
					}
					dtemp = c*c*d[i] + s*s*d[j] - 2.0*c*s*b;
					d[j] = s*s*d[i] + c*c*d[j] +  2.0*c*s*b;
					d[i] = dtemp;
				}  /* end if */
			} /* end for i */
		} /* end for j */
	} /* end for l */

	nrot = l;

	for (j = 0; j <= 2; j++) {
		k = j;
		dtemp = d[k];
		for (i = j+1; i <= 3; i++) {
			if(d[i] < dtemp) {
				k = i;
				dtemp = d[k];
			}
		}

		if(k > j) {
			d[k] = d[j];
			d[j] = dtemp;
			for (i = 0; i <= 3; i++) {
				dtemp = v[i][k];
				v[i][k] = v[i][j];
				v[i][j] = dtemp;
			}
		}
	}
}
