#pragma once

#ifndef _EXTERNAL_MATH_LIBS_PROTOTYPES_H
#define _EXTERNAL_MATH_LIBS_PROTOTYPES_H

extern "C"
{
	// LAPACK routines
	void dsptrd_(char*, int*,double*,double*,double*,double*,int*);
	void dstevr_(char*, char*, int*, double*, double*, double*, double*, int*, int*, double*, int*, double*,
						double*, int*, int* , double*, int*, int*, int*, int*);
	void dopmtr_(char*,char*,char*,int*,int*,double*,double*,double*,int*,double*,int*);
	void dsyev_(char*, char*,int*,double*,int*,double*,double*,int*,int*);

	// BLAS routines
	double dnrm2_(int *n, const double * const x, int *incx);
	double ddot_(int * n, const double * const x, int * incX, const double * const y, int * incY);
	void daxpy_(int * n, double * alpha, const double * const x, int * incx, double * const y, int * incy);
}


// Inlined wrappers for Blas and Lapack routines when using gcc

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the symetric sparse tridiagonal matrices with double precision.
///
/// \param char1 [In] U/L These values define if A matrix corresponds to Upper triangle or Lower triangle.
/// \param ln [In] Order of A matrix.
///
/// \param D [Out] Diagonal elements of the tridiagonal matrix computed. Dimension (ln).
/// \param E [Out] Super/Subdiagonal elements of the tridiagonal matrix computed. Dimension (ln-1).
/// \param tau [Out] The scalar factors of the elementary reflectors. Dimension (ln-1).
/// \param info [Out] 0= succesfull. <0 ilegal values.
///
/// \param A [In/Out] Upper/Lower triangle of the symetric matrix written in one dimension. Dimension (ln*(ln+1)/2).
///
/// \author icabeza
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline void dsptrd(char* char1, int* ln, double* A , double*D , double* E, double* tau, int* info)
{
	dsptrd_(char1, ln, A, D, E, tau, info);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the eigenvectors and eigenvalues in a range (rzer,ln)
/// of a tridiagonal matrix with double precision.
///
/// \param cha1 [In] Mode of calculations.
/// \param N [In] Compute eigenvalues.
/// \param V [In] Compute eigenvectors and eigenvalues.
/// \param cha2 [In] All eigenvalues will be found.
/// \param V [In] All eigenvalues in the half interval (rzer,rzer) will be found.
/// \param I [In] All eigenvalues between rzer and rzer will be found.
/// \param ln [In] Order of A matrix.
/// \param D [In] Diagonal elements of the tridiagonal matrix computed. Dimension (ln).
/// \param E [In] Super/Subdiagonal elements of the tridiagonal matrix computed. Dimension (ln-1).
/// \param rzer [In] Ignored in our choice (I).
/// \param rzer [In] Ignored in our choice (I).
/// \param eigs [In] Range to start to compute eigenvectors and eigenvalues.
/// \param eige [In] Range to end computation of eigenvectors and eigenvalues.
/// \param delta [In] Maximum tolerance allowed for the calculations of eigenvalues/eigenvectors.
/// \param lwork [In] Dimension of work.
/// \param liwork [In] Dimension of iwork.
///
/// \param nume [Out] Total number of eigenvalues found.
/// \param eig [Out] Array of eigenvalues found in ascending order.
/// \param eigv [Out] Matrix array of eigenvectors. Dimension (ln,max(1,nume))
/// \param ln [Out] Leading dimension of array eigv ---> Max(1,ln).
/// \param work [Out] Workspace/output. Dimension max(1,lwork).
/// \param iwork [Out] Workspace/output. Dimension max(1,liwork).
/// \param info [Out] 0= succesfull. <0 ilegal values.
///
/// \author icabeza
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline void dstevr(char* cha1, char* cha2, int * ln1, double* D, double* E, double * rzer1,
		double * rzer2, int* eigs, int* eige, double* delta, int* nume, double* eig,
		double* eigv, int * ln2, int* isuppz, double* work, int* lwork, int* iwork, int* liwork, int* info)
{
	dstevr_(cha1, cha2, ln1, D, E, rzer1, rzer2, eigs, eige, delta, nume, eig, eigv, ln2, isuppz, work, lwork, iwork, liwork, info);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the eigenvectors and eigenvalues of a matrix from eigenvectors and eigenvalues
/// of the tridiagonal matrix.
///
/// http://www.netlib.org/lapack/lapack-3.1.1/html/dopmtr.f.html
///
/// \param cha1 [In] These values define the operation form.
/// \param 'L' = Left, 'R' = Right
/// \param cha2 [In] Define if matrix is storage Upper or Lower
/// 'U' = Upper, 'L' = Lower
/// \param cha3 [In] Define if the matrix is transposed or not.
/// 'N' = Not transposed, 'T' = transposed
/// \param ln [In] Order of A matrix.
/// \param nume [In] Total number of eigenvalues.
/// \param A [In] Upper/Lower triangle of the symetric matrix written in one dimension. Dimension (ln*(ln+1)/2).
/// \param tau [In] The scalar factors of the elementary reflectors. Dimension (ln-1).
///
/// \param work [Out] Workspace/output. Dimension max(1,lwork).
/// \param info [Out] 0= succesfull. <0 ilegal values.
///
/// \param eigv [In/Out] Matrix array of eigenvectors. Dimension (ln,max(1,nume))
///
/// \author icabeza
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline void dopmtr(char* cha1, char* cha2,char* cha3, int* ln1, int* nume,
		double* A, double* tau, double* eigv, int* ln2, double* work, int* info)
{
	dopmtr_(cha1, cha2, cha3, ln1, nume, A, tau, eigv, ln2, work, info);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the eigenpairs of dense symetric matrices with double precision.
///
/// \param jobz [In] N/V Computes eigenvalues or eigenpairs, respectively
/// \param uplo [In] U/L The values defined in A matrix corresponds to Upper triangle or Lower triangle.
/// \param n [In] Order of A matrix.
/// \param lda [In] dimension of the array.
///
/// \param w [Out] Eigenvalues vector
/// \param work [Out] Temporal variable to work
/// \param lwork [Out] Size of work
/// \param info [Out] 0= succesfull. <0 ilegal values.
///
/// \param A [In/Out] Upper/Lower triangle of the symetric matrix written in one dimension. Dimension (ln*(ln+1)/2).
///
/// \author icabeza
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline void dsyev(char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info)
{
	dsyev_(jobz, uplo, n, a, lda, w, work, lwork, info);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the dot product of two vectors
///
/// \param n [In] Number of elements in each vector.
/// \param x [In] Array of dimension (n-1) * |incx| + 1.
/// Array x contains the first vector operand.
/// \param incx [In] Increment between elements of x.
/// If incx = 0, the results will be unpredictable.
/// \param y [In] Array of dimension (n-1) * |incy| + 1.
/// Array y contains the second vector operand.
/// \param incy [In] Increment between elements of y.  If incy = 0, the results will be unpredictable.
///
/// \return The dot product of the two vectors
///
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline double ddot(int * n, const double * const x, int * incX, const double * const y, int * incY)
{
	return ddot_(n, x, incX, y, incY);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes constant times a vector plus a vector: y  <--  alpha*x + y
///
/// \param n [In] Number  of  elements in the vectors.
/// If n <= 0, these routines return without any computation.
/// \param alpha [In] If alpha = 0 this routine returns without any computation.
/// \param x [In] Array of dimension (n-1) * |incx| + 1.
/// Contains the vector  to be scaled before summation.
/// \param incx [In] Increment between elements of x.
///	If incx = 0, the results will be unpredictable.
/// \param incy [In] Increment between elements of y.
/// If incy = 0, the results will be unpredictable.
///
/// \param y [In/Out] Array of dimension (n-1) * |incy| + 1.
///	Before calling the routine, y contains the vector to be summed.
///	After the routine ends, y contains the result of the summation.
///
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline void daxpy(int * n, double * sa, const double * const sx, int * incx, double * const sy, int * incy)
{
	daxpy_(n, sa, sx, incx, sy, incy);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function computes the euclidean norm of a vector
///
/// \param n [In] Number of elements of the vector
/// \param x [In] A one-dimensional array x of length at least (1+(n-1)*|incx|),
///				   containing the elements of the vector x.
/// \param incx [In] The increment for the array X.
///		If incx > 0, vector x is stored forward in the array, so that x(i) is stored in
///		location x(1+(i-1)*incx).
///		If incx < 0, vector x is stored backward in the array, so that x(i) is stored in
///		location x(1+(n-i)*|incx|).
///		If incx = 0, only the first element is accessed.
///
/// \return The euclidean norm of the vector
///
/// \author mrivero
/// \date 08/02/2011
///////////////////////////////////////////////////////////////
inline double dnrm2(int *n, const double * const x, int *incx)
{
	return dnrm2_(n, x, incx);
}


#endif /*#define _EXTERNAL_MATH_LIBS_PROTOTYPES_H*/
