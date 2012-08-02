
#include "RMSDTools.h"
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;


void RMSDTools::superpose(unsigned int n, double * const coord_fit, double* const coord_ref)
{
	double center_fit[3];
	superposeMatrix(n, coord_fit, coord_ref, center_fit);

	// Shift to the origin
	RMSDTools::shift3D(n, coord_fit, center_fit, +1.);
}


void RMSDTools::superposeMatrix(unsigned int n, double * const coord_fit, double* const coord_ref, double * const center_fit)
{
	double center_ref[3]={0,0,0};
	double u[3][3];
	double q[4];

	// Initialize center_fit as (0.0, 0.0, 0.0)
	center_fit[0] = 0.0;
	center_fit[1] = 0.0;
	center_fit[2] = 0.0;

	// Fit specified atoms of fit_molecule to those of ref_molecule

	// Get rotation matrix
	getRotationMatrixGetCentersAndShiftMolecules(center_fit, center_ref, n, coord_fit, coord_ref, u, q);

	// Rotate the reference molecule by the rotation matrix u
	RMSDTools::rotate3D(n, coord_ref, u);

	// Shift the reference molecule to the geometric center of the fit
	RMSDTools::shift3D(n, coord_ref, center_fit, +1.);
}

void RMSDTools::getRotationMatrixGetCentersAndShiftMolecules(double * const center_fit, double * const center_ref, unsigned int num_atoms, double * const coord_fit, double * const coord_ref, double u[][3], double * const q)
{
	// Geometrical center of reference coordinates
	RMSDTools::geometricCenter(num_atoms, coord_ref, center_ref);

	// Geometric center of fitted coordinates
	RMSDTools::geometricCenter(num_atoms, coord_fit, center_fit);

	// Shift reference coordinates to origin
	RMSDTools::shift3D(num_atoms, coord_ref, center_ref, -1.);

	// Shift fit coordinates to the origin
	RMSDTools::shift3D(num_atoms, coord_fit, center_fit, -1.);

	// Get rotation matrix
	RMSDTools::superpositionQuatFit(num_atoms, coord_ref, coord_fit, q, u);
}


void RMSDTools::geometricCenter(unsigned int n, const double * const x, double * const center)
{
	unsigned int i;

	// Initialize variables before the loop
	center[0] = 0.0;
	center[1] = 0.0;
	center[2] = 0.0;

	// Computed the weighted geometric center
	for(i=0; i<n; i++)
	{
		RMSDTools::ourDaxpy(3, 1., (x+3*i), center);
	}

	RMSDTools::multiplyVectorByScalar(center, 1.0/n, 3);
}

void RMSDTools::multiplyVectorByScalar(double * const x, double scalar, unsigned int n)
{
	for(unsigned int i=0; i<n; ++i)
	{
		*(x+i) *= scalar;
	}
}

void RMSDTools::multiplyVectorByScalar(double * y, const double * const x, double scalar, unsigned int n)
{
	for(unsigned int i=0; i<n; ++i)
	{
		*(y+i) = *(x+i) * scalar;
	}
}


void RMSDTools::ourDaxpy(int num_elems, double sa, const double * const sx, double * const sy)
{
	int inc = 1;
	daxpy(&num_elems, &sa, sx, &inc, sy, &inc);
}

void RMSDTools::shift3D(unsigned int numberOfPoints, double * const x, double trans[3], double scalar)
{
	const unsigned int num_coords = 3 * numberOfPoints;

	double shiftVector[3];
	RMSDTools::multiplyVectorByScalar(shiftVector, trans, scalar, 3);

	for(unsigned int i=0; i<num_coords; i+=3)
	{
		RMSDTools::shiftOne3dPoint(x+i, shiftVector);
	}
}

void RMSDTools::shiftOne3dPoint(double * const x, double trans[3])
{
	for(unsigned int i=0; i<3; ++i)
	{
		*(x+i) += trans[i];
	}
}

void RMSDTools::superpositionQuatFit(unsigned int n, const double * const x, const double * const y, double q[4], double u[3][3])
{
	vector<double> w(n, 1.0);

	superpositionQuatFit(n, x, y, &w[0], q, u);
}

void RMSDTools::superpositionQuatFit(unsigned int n, const double * const x, const double * const y, const double * const w, double q[4], double u[3][3])
{
	vector<double> upperQuadMatrix(16, 0);

	// Generate the upper triangle of the quadratic matrix
	generateUpperQuadraticMatrix(&upperQuadMatrix[0], n, y, x, w);

	// Compute quaternion
	computeFittedQuaternionFromUpperQuadraticMatrix(&upperQuadMatrix[0], q);

	// Generate the rotation matrix
	generateLeftRotationMatrixFromNormalizedQuaternion(q, u);
}

void RMSDTools::generateUpperQuadraticMatrix(double * const upperQuadMatrix, unsigned int n, const double * const y, const double *const x, const double * const w)
{
	//Initialize upperQuadMatrix
	for(unsigned int i=0; i<16; ++i)
	{
		upperQuadMatrix[i] = 0.0;
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

		xxyx += x[k] * y[k] * w[i];
		xxyy += x[k] * y[k+1] * w[i];
		xxyz += x[k] * y[k+2] * w[i];
		xyyx += x[k+1] * y[k] * w[i];
		xyyy += x[k+1] * y[k+1] * w[i];
		xyyz += x[k+1] * y[k+2] * w[i];
		xzyx += x[k+2] * y[k] * w[i];
		xzyy += x[k+2] * y[k+1] * w[i];
		xzyz += x[k+2] * y[k+2] * w[i];
	}

	upperQuadMatrix[0] = xxyx + xyyy + xzyz;

	upperQuadMatrix[4] = xzyy - xyyz;
	upperQuadMatrix[5] = xxyx - xyyy - xzyz;

	upperQuadMatrix[8] = xxyz - xzyx;
	upperQuadMatrix[9] = xxyy + xyyx;
	upperQuadMatrix[10] = xyyy - xzyz - xxyx;

	upperQuadMatrix[12] = xyyx - xxyy;
	upperQuadMatrix[13] = xzyx + xxyz;
	upperQuadMatrix[14] = xyyz + xzyy;
	upperQuadMatrix[15] = xzyz - xxyx - xyyy;
}

void RMSDTools::computeFittedQuaternionFromUpperQuadraticMatrix(double * const upperQuadMatrix, double q[4])
{
	// Diagonalize upperQuadMatrix
	int info;
	int dimension = 4;
	int lda = 4;

	vector<double> work1(4, 0.0);

	// 136 is the optimal size for LWORK
	int lwork = 136;
	vector<double> work2(lwork, 0.0);

	// Get eigenvalues and eigenVectors
	ourDsyevUpperVector(&dimension, upperQuadMatrix, &lda, &work1[0], &work2[0], &lwork, &info);

	// Resulting quaternion
	q[0] = upperQuadMatrix[12];
	q[1] = upperQuadMatrix[13];
	q[2] = upperQuadMatrix[14];
	q[3] = upperQuadMatrix[15];
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

void RMSDTools::ourDsyevUpperVector(int *n, double * eigenVectors, int *lda, double * eigenValues,
									double * work, int * lwork, int * info)
{
	char vector='V';
	char upper='U';
	dsyev(&vector, &upper, n, eigenVectors, lda, eigenValues, work, lwork, info);
}

void RMSDTools::rotate3D(unsigned int numberOfPoints, const double * const x, double * const y, double u[3][3])
{
	const unsigned int num_coords = 3 * numberOfPoints;

	// We go through all points
	for(unsigned int i=0; i<num_coords; i+=3)
	{
		// An rotate each of them
		RMSDTools::rotate3DPoint( (x+i), (y+i), u);
	}
}

void RMSDTools::rotate3D(unsigned int n, double * const x, double u[3][3])
{
	const unsigned int num_coords = 3 * n;

	double * y = new double[num_coords];

	rotate3D(n, x, y, u);

	std::copy(y, y + num_coords, x);

	if (y!= NULL)
	delete [] y;
}

void RMSDTools::rotate3DPoint(const double * const x, double * const y, double u[3][3])
{
	*y = u[0][0] * *x + u[0][1] * *(x+1) + u[0][2] * *(x+2);
	*(y+1) = u[1][0] * *x + u[1][1] * *(x+1) + u[1][2] * *(x+2);
	*(y+2) = u[2][0] * *x + u[2][1] * *(x+1) + u[2][2] * *(x+2);
}

double RMSDTools::calcRMS(const double * const x, const double * const y, unsigned int num_atoms)
{
	double sum_res = 0.0;

	for(unsigned int i=0; i<num_atoms*3; ++i)
	{
		sum_res += (x[i] - y[i]) * (x[i] - y[i]);
	}

	return sqrt(sum_res/num_atoms);
}
