#ifndef _MATH_TOOLS_H
#define _MATH_TOOLS_H

#include "ExternalMathLibsPrototypes.h"

namespace RMSDTools{

	void superpose(unsigned int n, double * const coord_fit, double* const coord_ref);
	
	void superposeMatrix(unsigned int n, double * const coord_fit, double* const coord_ref, double * const center_fit);
	
	void getRotationMatrixGetCentersAndShiftMolecules(double * const center_fit, double * const center_ref, unsigned int num_atoms, double * const coord_fit, double * const coord_ref, double u[][3], double * const q);
	
	void geometricCenter(unsigned int n, const double * const x, double * const center);
	
	void multiplyVectorByScalar(double * const x, double scalar, unsigned int n);
	
	void multiplyVectorByScalar(double * y, const double * const x, double scalar, unsigned int n);
	
	void ourDaxpy(int num_elems, double sa, const double * const sx, double * const sy);
	
	void shift3D(unsigned int numberOfPoints, double * const x, double trans[3], double scalar);
	
	void shiftOne3dPoint(double * const x, double trans[3]);
	
	void superpositionQuatFit(unsigned int n, const double * const x, const double * const y, double q[4], double u[3][3]);
	
	void superpositionQuatFit(unsigned int n, const double * const x, const double * const y, const double * const w, double q[4], double u[3][3]);
	
	void generateUpperQuadraticMatrix(double * const upperQuadMatrix, unsigned int n, const double * const y, const double *const x, const double * const w);
	
	void computeFittedQuaternionFromUpperQuadraticMatrix(double * const upperQuadMatrix, double q[4]);
	
	void generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3]);
	
	void ourDsyevUpperVector(int *n, double * eigenVectors, int *lda, double * eigenValues,
										double * work, int * lwork, int * info);
										
	void rotate3D(unsigned int numberOfPoints, const double * const x, double * const y, double u[3][3]);

	void rotate3D(unsigned int n, double * const x, double u[3][3]);
	
	void rotate3DPoint(const double * const x, double * const y, double u[3][3]);
	
	double calcRMS(const double * const x, const double * const y, unsigned int num_atoms);
}

#endif