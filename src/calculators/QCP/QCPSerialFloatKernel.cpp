#include <cmath>
#include <iostream>
#include "QCPSerialFloatKernel.h"
#include "../RMSDTools.h"
#include "../RMSDCalculationData.h"

using namespace std;
/**
 * Highly underperforming class to evaluate CUDA float version.
 */


///////////////////////////////////////////////////////////////
/// Implements the 'inner product' operation of Douglas Theobald QCP superposition method (see : http://theobald.brandeis.edu/qcp/
/// and "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial."  Acta Crystallogr A 61(4):478-480
/// for more info).
///
/// \param 	A [In/Out] A 3x3 matrix (                blank space for a real explanation of this parameter)
///
/// \param 	first_conformation_coords [In] Array containing the coordinates of the reference conformation.
///
/// \param 	second_conformation_coords [In] Array containing the coordinates of the conformation to be measured.
///
/// \param 	number_of_atoms [In] Number of atoms of both conformations (must be the same, of course).
///
/// \return The E0 parameter (upper bound for max Eigenvalue).
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
float QCPSerialFloatKernel::innerProduct(float* A,
												float* first_conformation_coords,
												float* second_conformation_coords,
												int number_of_atoms){
    float x1, x2, y1, y2, z1, z2;
    int i;
    float G1 = 0.0, G2 = 0.0;

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0;

	int total_number_of_coordinates = 3* number_of_atoms;

    for (i = 0; i < total_number_of_coordinates; i+=3){
        x1 = first_conformation_coords[i];
        y1 = first_conformation_coords[i+1];
        z1 = first_conformation_coords[i+2];

        x2 = second_conformation_coords[i];
        y2 = second_conformation_coords[i+1];
        z2 = second_conformation_coords[i+2];

        G1 += x1 * x1 + y1 * y1 + z1 * z1;

        G2 += x2 * x2 + y2 * y2 + z2 * z2;

        A[0] +=  (x1 * x2);
        A[1] +=  (x1 * y2);
        A[2] +=  (x1 * z2);

        A[3] +=  (y1 * x2);
        A[4] +=  (y1 * y2);
        A[5] +=  (y1 * z2);

        A[6] +=  (z1 * x2);
        A[7] +=  (z1 * y2);
        A[8] +=  (z1 * z2);
    }

    return (G1 + G2) * 0.5;
}

///////////////////////////////////////////////////////////////
/// \remarks
///	This function ports the second part of Douglas Theobald QCP superposition method (see 'innerProduct').
///
/// \param 	A [In] A 3x3 matrix (more ?).
///
/// \param 	E0 [In] Upper bound for the maximum eigenvalue (?).
///
/// \param 	number_of_atoms [In] Number of atoms of conformations used to get A and E0.
///
/// \return The actual rmsd between both conformations.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
float QCPSerialFloatKernel::calcRMSDForTwoConformationsWithTheobaldMethod(
		float *A,
		float E0,
		int number_of_atoms,
		float* rot_matrix){
    float Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    float Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    float C[4];
    float mxEigenV;
    float b, a, delta, x2;
    float oldg = 0.0;
    float evalprec = 1e-11;

    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    mxEigenV = E0;
    for (int i = 0; i < 50; ++i){
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
            break;
    }

    if (rot_matrix != NULL ){
		float a11, a12, a13, a14, a21, a22, a23, a24,
				a31, a32, a33, a34, a41, a42, a43, a44;

		float a3344_4334, a3244_4234, a3243_4233,
				a3143_4133, a3144_4134, a3142_4132;

		float q1, q2, q3, q4, normq;

		float evecprec = 1e-6;

		float a2, x2, y2, z2;

		float xy, az, zx, ay, yz, ax;

		a11 = SxxpSyy + Szz - mxEigenV;
		a12 = SyzmSzy;
		a13 = -SxzmSzx;
		a14 = SxymSyx;
		a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
		a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
		a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
		a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
		a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
		a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
		q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
		q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
		q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
		q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

		float qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

		if (qsqr < evecprec)
		{
			q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
			q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
			q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
			q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
			qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

			if (qsqr < evecprec)
			{
				float a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
				float a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
				float a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

				q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
				q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
				q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
				q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
				qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

				if (qsqr < evecprec)
				{
					q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
					q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
					q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
					q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
					qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;

					/*if (qsqr < evecprec)
					{
						// if qsqr is still too small, return the identity matrix.
						rot_matrix[0] = rot_matrix[4] = rot_matrix[8] = 1.0;
						rot_matrix[1] = rot_matrix[2] = rot_matrix[3] = rot_matrix[5] = rot_matrix[6] = rot_matrix[7] = 0.0;
						//cout<<"Qsqr is too small!"<<endl;
					}*/
				}
			}
		}

		normq = sqrt(qsqr);
		q1 /= normq;
		q2 /= normq;
		q3 /= normq;
		q4 /= normq;

		a2 = q1 * q1;
		x2 = q2 * q2;
		y2 = q3 * q3;
		z2 = q4 * q4;

		xy = q2 * q3;
		az = q1 * q4;
		zx = q4 * q2;
		ay = q1 * q3;
		yz = q3 * q4;
		ax = q1 * q2;

		rot_matrix[0] = a2 + x2 - y2 - z2;
		rot_matrix[1] = 2 * (xy + az);
		rot_matrix[2] = 2 * (zx - ay);
		rot_matrix[3] = 2 * (xy - az);
		rot_matrix[4] = a2 - x2 + y2 - z2;
		rot_matrix[5] = 2 * (yz + ax);
		rot_matrix[6] = 2 * (zx + ay);
		rot_matrix[7] = 2 * (yz - ax);
		rot_matrix[8] = a2 - x2 - y2 + z2;
    }

    if (isnan(mxEigenV)){
    	return 0.0;
    }
    else{
    	return sqrt(fabs(2.0 * (E0 - mxEigenV)/number_of_atoms));
    }
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Wrapping function for Douglas Theobald QCP superposition method to calculate the RMSD for two conformations.
///
/// \param 	first_conformation_coords [In] Array containing the coordinates of the reference conformation.
///
/// \param 	second_conformation_coords [In] Array containing the coordinates of the conformation to be measured.
///
/// \param 	number_of_atoms [In] Number of atoms of both conformations (must be the same, of course).
///
/// \return The actual rmsd between both conformations.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
float QCPSerialFloatKernel::calcRMSDOfTwoConformations(
		float* first_conformation_coords,
		float* second_conformation_coords,
		int number_of_atoms,
		float* rot_matrix){

	float A[9];
	float E0 = innerProduct(A, first_conformation_coords, second_conformation_coords, number_of_atoms);
	return calcRMSDForTwoConformationsWithTheobaldMethod(A, E0, number_of_atoms, rot_matrix);
}

///////////////////////////////////////////////////////////////
/// Given a trajectory with N conformations stored in 'all_coordinates', it calculates the rmsd between the conformation
/// with 'i' = 'base_conformation_id' and the conformations in the range [j+1,M)
///
/// \param 	all_coordinates [In] The array storing all coordinates for all conformations in order:
/// (conf0_atom_0_x,conf0_atom_0_y,conf0_atom_0_z,conf0_atom_1_x,conf0_atom_1_y ... confN-1_atom_N_z)
///
/// \param 	base_conformation_id [In] The reference conformation id ('i' in the above explanation).
///
/// \param 	other_conformations_starting_id [In] The value of j in the above explanation.
///
/// \param 	number_of_conformations [In] Number of conformations stored in the coordinates array.
///
/// \param 	number_of_atoms [In] Number of atoms of any of the conformations.
///
/// \param 	rmsd [In] Array where the rmsds are stored (rmsd[j]  is the rmsd between conformation 'i' and
///	conformation 'j').
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
void QCPSerialFloatKernel::oneVsFollowingFitEqualCalcCoords(
		double* reference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	float rot_matrix[9];

	float* reference_tmp = new float[data->fittingConformationLength];
	float* second_conf_tmp = new float[data->fittingConformationLength];
	for(int i =0; i < data->fittingConformationLength;++i){
		reference_tmp[i] = (float) reference[i];
	}

	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double* second_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		for(int i =0; i < data->fittingConformationLength;++i){
			second_conf_tmp[i] = (float) second_conformation_coords[i];
		}

		RMSDTools::initializeTo(rot_matrix, 0.0, 9);

		float rmsd_val = calcRMSDOfTwoConformations(	reference_tmp,
														second_conf_tmp,
														data->atomsPerFittingConformation,
														rot_matrix);
		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = (double) rmsd_val;
		}

		RMSDTools::rotate3D(data->atomsPerFittingConformation, second_conf_tmp, rot_matrix);
		for(int i =0; i < data->fittingConformationLength;++i){
			second_conformation_coords[i] = (float) second_conf_tmp[i];
		}
	}

	delete [] reference_tmp;
	delete [] second_conf_tmp;
}

void QCPSerialFloatKernel::oneVsFollowingFitDiffersCalcCoords(
		double* fitReference,
		double* calcReference,
		int reference_conformation_number,
		double* rmsd,
		RMSDCalculationData* data){

	float rot_matrix[9];

	float* fit_reference_tmp = new float[data->fittingConformationLength];
	float* calc_reference_tmp = new float[data->calculationConformationLength];
	for(int i =0; i < data->fittingConformationLength;++i){
		fit_reference_tmp[i] = (float) fitReference[i];
	}
	if (calcReference!=NULL){
		for(int i =0; i < data->calculationConformationLength;++i){
			calc_reference_tmp[i] = (float) calcReference[i];
		}
	}

	float* fit_second_conf_tmp = new float[data->fittingConformationLength];
	float* calc_second_conf_tmp = new float[data->calculationConformationLength];

	for (int second_conformation_index = reference_conformation_number + 1;
			second_conformation_index < data->numberOfConformations; ++second_conformation_index){

		double* fit_conformation_coords = data->getFittingConformationAt(second_conformation_index);
		for(int i =0; i < data->fittingConformationLength;++i){
			fit_second_conf_tmp[i] = (float) fit_conformation_coords[i];
		}

		double* calc_conformation_coords =  data->getCalculationConformationAt(second_conformation_index);
		for(int i =0; i < data->calculationConformationLength;++i){
			calc_second_conf_tmp[i] = (float) calc_conformation_coords[i];
		}

		RMSDTools::initializeTo(rot_matrix, 0.0, 9);

		calcRMSDOfTwoConformations(	fit_reference_tmp,
									fit_second_conf_tmp,
									data->atomsPerFittingConformation,
									rot_matrix);

		RMSDTools::rotate3D(data->atomsPerFittingConformation,
				fit_second_conf_tmp,
				rot_matrix);
		RMSDTools::rotate3D(data->atomsPerCalculationConformation,
				calc_second_conf_tmp,
				rot_matrix);

		if(rmsd!=NULL){
			rmsd[second_conformation_index-(reference_conformation_number+1)] = RMSDTools::calcRMS(
					calc_reference_tmp,
					calc_second_conf_tmp,
					data->atomsPerCalculationConformation);
		}

		for(int i =0; i < data->fittingConformationLength;++i){
			fit_conformation_coords[i] = (float) fit_second_conf_tmp[i];
		}

		for(int i =0; i < data->calculationConformationLength;++i){
			calc_conformation_coords[i] = (float) calc_second_conf_tmp[i];
		}
	}

	delete [] fit_reference_tmp;
	delete [] calc_reference_tmp;
	delete [] fit_second_conf_tmp;
	delete [] calc_second_conf_tmp;

}

