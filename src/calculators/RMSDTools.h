#ifndef _MATH_TOOLS_H
#define _MATH_TOOLS_H

#include <vector>
#include <utility>
#include "symmGroups.h"

namespace RMSDTools{

	void superpose(unsigned int n, double * const coord_fit, double* const coord_ref);

	void superpose(unsigned int fit_n, double * const fit_coords, double* const fit_ref_coords,
			unsigned int rmsd_n, double* const calc_coords);
	
	void centerAllAtOrigin(unsigned int atomsPerConformation, unsigned int numberOfConformations,
					double * const all_coords, double* const translations);

	void centerAllAtOrigin(unsigned int atomsPerConformation, unsigned int numberOfConformations,
						double * const all_coords);

	void applyTranslationsToAll(unsigned int atomsPerConformation, unsigned int numberOfConformations,
			double * const all_coords, double* const translations, int sign = 1);

	void applyTranslationToAll(unsigned int atomsPerConformation, unsigned int numberOfConformations,
			double * const all_coords, double* const translation_vector);

	void geometricCenter(unsigned int n, const double * const x, double * const center);
	
	void translate(unsigned int numberOfPoints, double * const x, double trans[3], double scalar);
	
	void superpositionQuatFit(unsigned int n, const double * const x, const double * const y, double q[4], double u[3][3]);
	
	void generateUpperQuadraticMatrix(double upperQuadMatrix[4][4], unsigned int n, const double * const y, const double *const x);
	
	void computeFittedQuaternionFromUpperQuadraticMatrix(double upperQuadMatrix[4][4], double q[4]);
	
	void generateLeftRotationMatrixFromNormalizedQuaternion(double q[4], double u[3][3]);
	
	void rotate3D(unsigned int n, double * const x, double*);

	void rotate3D(unsigned int n, float * const x, float*);

	void rotate3D(unsigned int n, double * const x, double u[3][3]);
	
	void rotate3D(unsigned int n, float * const x, float u[3][3]);

	double calcRMS(const double * const x, const double * const y, unsigned int num_atoms);

	double calcRMS(const float * const x, const float * const y, unsigned int num_atoms);

	void jacobi(double a[4][4], double d[4], double v[4][4], int nrot = 30);

	void transposeMatrix(const double (*const source_matrix)[3], double matrix[3][3]);

	void transposeMatrix(const double (*const source_matrix)[4], double matrix[3][3]);

	bool diagonalize_symmetric(double matrix[3][3], double eigen_vec[3][3], double eigenval[3]);

	void initializeTo(double*, double, int);

	void initializeTo(float*, float, int);

	void copyArrays(double*, double*, int);

	void calculateMeanCoordinates(double* , double* , int , int );

	void normalize(double*);

	double dot(double*, double*);

	void cross(double* , double*, double*);

	bool jacobi3(double a[3][3], double d[3], double v[3][3], int n_rot = 50);

	void swap_atoms(double* coordinates, int atom_i, int atom_j);

	void applySymmetryGroup(double* coordinates, std::vector<std::pair<int,int> >& symm_group);

	void calcRecursiveSymmGroupApplication(double* reference, double* superposed_conformation,
												int number_of_atoms, symmGroups* symm_groups,
												int applied_symm_group, std::vector<double>& rmsds);

	double calcMinRMSDOfAllSymmetryGroups(	double* reference, double* superposed_conformation,
												int number_of_atoms, symmGroups* symm_groups);

}

#endif
