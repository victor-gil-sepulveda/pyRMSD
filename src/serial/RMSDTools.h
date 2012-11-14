#ifndef _MATH_TOOLS_H
#define _MATH_TOOLS_H

namespace RMSDTools{

	void superpose(unsigned int n, double * const coord_fit, double* const coord_ref);
	
	void superposeMatrix(unsigned int n, double * const coord_fit, double* const coord_ref, double * const center_fit);
	
	void getRotationMatrixGetCentersAndShiftMolecules(double * const center_fit, double * const center_ref, unsigned int num_atoms, double * const coord_fit, double * const coord_ref, double u[][3], double * const q);
	
	void geometricCenter(unsigned int n, const double * const x, double * const center);
	
	void shift3D(unsigned int numberOfPoints, double * const x, double trans[3], double scalar);
	
	void superpositionQuatFit(unsigned int n, const double * const x, const double * const y, double q[4], double u[3][3]);
	
	void generateUpperQuadraticMatrix(double upperQuadMatrix[4][4], unsigned int n, const double * const y, const double *const x);
	
	void computeFittedQuaternionFromUpperQuadraticMatrix(double upperQuadMatrix[4][4], double q[4]);
	
	void generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3]);
	
	void rotate3D(unsigned int n, double * const x, double u[3][3]);
	
	double calcRMS(const double * const x, const double * const y, unsigned int num_atoms);

	void jacobi(double a[4][4], double d[4], double v[4][4], int nrot = 30);
}

#endif
