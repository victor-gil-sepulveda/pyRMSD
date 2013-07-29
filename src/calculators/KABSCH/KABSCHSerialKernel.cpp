/*
 * KABSCHSerialKernel.cpp
 *
 *  Created on: 07/03/2013
 *      Author: victor
 */

#include "KABSCHSerialKernel.h"
#include "../RMSDTools.h"
#include <cmath>
#include <iostream>
#include "../RMSDCalculationData.h"
using namespace std;

KABSCHSerialKernel::KABSCHSerialKernel() {}

KABSCHSerialKernel::~KABSCHSerialKernel() {}

void KABSCHSerialKernel::oneVsFollowingFitEqualCalcCoords(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			RMSDCalculationData* data){

		double rot_matrix[3][3];

		for (int second_conformation_index = reference_conformation_number + 1;
				second_conformation_index < data->numberOfConformations; ++second_conformation_index){

			double* second_conformation_coords = data->getFittingConformationAt(second_conformation_index);

			RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

			double rmsd_val = this->calculate_rotation_rmsd(
					reference,
					second_conformation_coords,
					data->atomsPerFittingConformation,
					rot_matrix);

			if(rmsd!=NULL){
				rmsd[second_conformation_index-(reference_conformation_number+1)] = rmsd_val;
			}

			RMSDTools::rotate3D(data->atomsPerFittingConformation,
					second_conformation_coords,
					rot_matrix);
		}
}

void KABSCHSerialKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	double rot_matrix[3][3];

	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double* fit_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		double* calc_conformation_coords =  data->getCalculationConformationAt(second_conformation_index);

		RMSDTools::initializeTo(rot_matrix[0], 0.0, 9);

		this->calculate_rotation_rmsd(
						fitReference,
						fit_conformation_coords,
						data->atomsPerFittingConformation,
						rot_matrix);

		RMSDTools::rotate3D(data->atomsPerFittingConformation,
				fit_conformation_coords,
				rot_matrix);
		RMSDTools::rotate3D(data->atomsPerCalculationConformation,
				calc_conformation_coords,
				rot_matrix);

		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(calcReference,
					calc_conformation_coords,
					data->atomsPerCalculationConformation);
		}
	}
}



//-------------------------------------------
// Code based on the work of Dr. Bosco K. Ho
//
// http://boscoh.com/code/rmsd.c
//
//-------------------------------------------
void KABSCHSerialKernel::calc_correlation_matrix_and_E0(
		double (*const R)[3],
		double* const _E0,
		const double* const referenceCoords,
		const double* const fitCoords,
		int numberOfAtoms){

	// Init R to 0
	double E0 = 0.0;

	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			R[i][j] = 0.0;
		}
	}

	for (int n = 0; n < numberOfAtoms; ++n){
		int coords_offset = n*3;
		/*
		* E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n)
		*/
		E0 +=  	fitCoords[coords_offset] * fitCoords[coords_offset] + referenceCoords[coords_offset] * referenceCoords[coords_offset] +
				fitCoords[coords_offset+1] * fitCoords[coords_offset+1] + referenceCoords[coords_offset+1] * referenceCoords[coords_offset+1] +
				fitCoords[coords_offset+2] * fitCoords[coords_offset+2] + referenceCoords[coords_offset+2] * referenceCoords[coords_offset+2];

		/*
		* correlation matrix R:
		*   R[i,j) = sum(over n): y(n,i) * x(n,j)
		*   where x(n) and y(n) are two vector sets
		*/
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				R[i][j] += fitCoords[coords_offset+i] * referenceCoords[coords_offset+j];
			}
		}
	}
	E0 *= 0.5;
	*_E0 = E0;
}

/*
 * calculate_rotation_matrix()
 *
 *   calculates the rotation matrix U and the
 * rmsd from the R matrix and E0:
 */
bool KABSCHSerialKernel::calculate_rotation_matrix(
		const double (*const R)[3],
		double (*const U)[3],
		double E0,
		double* residual){

	int i, j, k;
	double Rt[3][3], RtR[3][3];
	double left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
	double v[3];
	double sigma;

	// build Rt, transpose of R
	RMSDTools::transposeMatrix(R,Rt);

	// make symmetric RtR = Rt X R
	for (i=0; i<3; i++){
		for (j=0; j<3; j++)	{
			RtR[i][j] = 0.0;
			for (k = 0; k<3; k++){
				RtR[i][j] += Rt[k][i] * R[j][k];
			}
		}
	}

	if (!RMSDTools::diagonalize_symmetric(RtR, right_eigenvec, eigenval)){
		return(false);
	}

	//cout<<"Reigenvec "<<right_eigenvec[0][0]<<endl;
	/* right_eigenvec's should be an orthogonal system but could be left
	 * or right-handed. Let's force into right-handed system.
	 */
	RMSDTools::cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);

	/* From the Kabsch algorithm, the eigenvec's of RtR
	 * are identical to the right_eigenvec's of R.
	 * This means that left_eigenvec = R x right_eigenvec
	 */
	for (i=0; i<3; i++){
		for (j=0; j<3; j++){
			left_eigenvec[i][j] = RMSDTools::dot(&right_eigenvec[i][0], &Rt[j][0]);
		}
	}

	for (i=0; i<3; i++){
		RMSDTools::normalize(&left_eigenvec[i][0]);
	}

	//cout<<"Leigenvec "<<left_eigenvec[0][0]<<endl;
	/*
	 * Force left_eigenvec[2] to be orthogonal to the other vectors.
	 * First check if the rotational matrices generated from the
	 * orthogonal eigenvectors are in a right-handed or left-handed
	 * co-ordinate system - given by sigma. Sigma is needed to
	 * resolve this ambiguity in calculating the RMSD.
	 */
	RMSDTools::cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);

	if(RMSDTools::dot(v, &left_eigenvec[2][0]) < 0.0){
		sigma = -1.0;
	}
	else{
		sigma = 1.0;
	}

	for (i=0; i<3; i++){
		left_eigenvec[2][i] = v[i];
	}

	/* calc optimal rotation matrix U that minimizes residual */
	for (i=0;i<3; i++){
		for (j=0; j<3; j++){
			U[i][j] = 0.0;
			for (k=0; k<3; k++){
				//cout<<"U["<<i<<"]["<<j<<"] "<<left_eigenvec[k][i] <<" "<< right_eigenvec[k][j]<<endl;
				if(!isnan(left_eigenvec[k][i]) && !isnan(right_eigenvec[k][i]))
					U[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
			}
			//cout<<"U["<<i<<"]["<<j<<"] "<<U[i][j]<<endl;
		}
	}

	*residual = E0 - (double) sqrt(fabs(eigenval[0]))
						 - (double) sqrt(fabs(eigenval[1]))
						 - sigma * (double) sqrt(fabs(eigenval[2]));

	return true;
}

double KABSCHSerialKernel::calculate_rotation_rmsd(
		const double* const referenceCoords,
		const double* const fitCoords,
		int numberOfAtoms,
		double (*const U)[3]){

	double E0, residual;
	double R[3][3];

	this->calc_correlation_matrix_and_E0(
		&(R[0]),
		&E0,
		referenceCoords,
		fitCoords,
		numberOfAtoms);

	this->calculate_rotation_matrix(R, U, E0, &residual);

	residual = fabs(residual); /* avoids the awkward case of -0.0 */

	return sqrt(fabs((double) (residual)*2.0/((double) numberOfAtoms)));
}
