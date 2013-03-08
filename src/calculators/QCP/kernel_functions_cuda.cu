#include "kernel_functions_cuda.h"
#include <cstddef>

///////////////////////////////////////////////////////////////
/// \remarks
/// Uses one thread to center the coordinates an arbitrary large number of conformations based on the GPU block topology chosen.
///
/// \param 	number_of_conformations [In] The number of conformations stored in all_coordinates.
///
/// \param 	number_of_atoms [In] The number of atoms PER CONFORMATION.
///
/// \param 	all_coordinates [In/Out] An array containing all the coordinates of all the conformations of the
///	trajectory, where the first conformation first atom coordinates are x = all_coordinates[0],
/// y = all_coordinates[1], z = all_coordinates[2]; Second conformation's first atom's coordinates would be
/// x = all_coordinates[number_of_atoms*3 + 0], y = all_coordinates[number_of_atoms*3 + 1] ... and so on.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
__global__ void centerCoordsOfAllConformations(const int number_of_conformations, const int number_of_atoms, floating_point_type* all_coordinates){
	int current_conformation_id = blockDim.x*blockIdx.x + threadIdx.x;
	while(current_conformation_id < number_of_conformations){
		int coordinates_per_conformation = number_of_atoms * 3;
		floating_point_type* conformation_coordinates = &(all_coordinates[current_conformation_id*coordinates_per_conformation]);
		centerCoords(conformation_coordinates, number_of_atoms);
		current_conformation_id += blockDim.x*gridDim.x;
	}
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Uses one thread to center the coordinates of ONE conformation.
///
/// \param 	conformation_coordinates [In] Pointer to the starting point of this conformation's coordinates inside
/// the array containing ALL the coordinates.
///
/// \param 	number_of_atoms [In] The number of atoms of this conformation.
///
/// \author victor_gil
/// \date 05/10/2012
///////////////////////////////////////////////////////////////
__device__ void centerCoords( floating_point_type* conformation_coordinates, const int number_of_atoms){
    floating_point_type xsum = 0;
    floating_point_type ysum = 0;
    floating_point_type zsum = 0;

	int total_number_of_coordinates = 3* number_of_atoms;
    for (int i = 0; i < total_number_of_coordinates; i+=3){
        xsum += conformation_coordinates[i];
        ysum += conformation_coordinates[i+1];
        zsum += conformation_coordinates[i+2];
    }

    for (int i = 0; i < total_number_of_coordinates; i+=3){
        conformation_coordinates[i] -= xsum / number_of_atoms;
        conformation_coordinates[i+1] -= ysum / number_of_atoms;
        conformation_coordinates[i+2] -= zsum / number_of_atoms;
    }
}


///////////////////////////////////////////////////////////////
/// \remarks
/// Implements the 'inner product' operation of Douglas Theobald QCP superposition method (see : http://theobald.brandeis.edu/qcp/
/// and "Rapid calculation of RMSDs using a quaternion-based characteristic polynomial."  Acta Crystallogr A 61(4):478-480
/// for more info ).
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
__device__ floating_point_type innerProduct(floating_point_type* A, 
											floating_point_type* first_conformation_coords, 
											floating_point_type* second_conformation_coords, 
											const int number_of_atoms){
    floating_point_type  x1, x2, y1, y2, z1, z2;
    floating_point_type  G1 = 0.0f, G2 = 0.0f;

    A[0] = A[1] = A[2] = A[3] = A[4] = A[5] = A[6] = A[7] = A[8] = 0.0f;

	int total_number_of_coordinates = 3* number_of_atoms;
    
    for (int i = 0; i < total_number_of_coordinates; i+=3){
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
__device__ floating_point_type calcRMSDForTwoConformationsWithTheobaldMethod(
		floating_point_type *A,
		const floating_point_type E0,
		const int number_of_atoms,
		floating_point_type * rot_matrix){
    floating_point_type Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    floating_point_type Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    floating_point_type C[4];
    floating_point_type mxEigenV; 
    floating_point_type b, a, delta, x2;
    floating_point_type oldg = 0.0f;
    floating_point_type evalprec = 1e-11f;

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

    SyzSzymSyySzz2 = 2.0f*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0f * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0f * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

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
        delta = ((a*mxEigenV + C[0])/(2.0f*x2*mxEigenV + b + a));
        mxEigenV -= delta;
        if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
            break;
    }

    if (rot_matrix != NULL ){
    		double a11, a12, a13, a14, a21, a22, a23, a24,
    				a31, a32, a33, a34, a41, a42, a43, a44;

    		double a3344_4334, a3244_4234, a3243_4233,
    				a3143_4133, a3144_4134, a3142_4132;

    		double q1, q2, q3, q4, normq;

    		double evecprec = 1e-6;

    		double a2, x2, y2, z2;

    		double xy, az, zx, ay, yz, ax;

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

    		double qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

    		if (qsqr < evecprec)
    		{
    			q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
    			q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
    			q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
    			q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
    			qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

    			if (qsqr < evecprec)
    			{
    				double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
    				double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
    				double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

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
    return sqrt(fabs(2.0f * (E0 - mxEigenV)/number_of_atoms));
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
__device__ floating_point_type calcRMSDOfTwoConformations(
		floating_point_type* first_conformation_coords,
		floating_point_type* second_conformation_coords,
		const int number_of_atoms,
		floating_point_type* rot_matrix){

	floating_point_type A[9];
	floating_point_type E0 = innerProduct(A, first_conformation_coords, second_conformation_coords, number_of_atoms);
	return calcRMSDForTwoConformationsWithTheobaldMethod(
			A,
			E0,
			number_of_atoms,
			rot_matrix);
}

__device__ void rotate3D(
		unsigned int number_of_atoms,
		floating_point_type * const coords,
		floating_point_type* u){

	// We go through all selected atoms
	for(unsigned int i=0; i<number_of_atoms; ++i){
		int offset = i*3;
		floating_point_type x_tmp_0,x_tmp_1,x_tmp_2;
		x_tmp_0 = coords[offset];
		x_tmp_1 = coords[offset+1];
		x_tmp_2 = coords[offset+2];

		// An rotate each of them
		coords[offset] 	= u[0] * x_tmp_0 + u[1] * x_tmp_1 + u[2] * x_tmp_2;
		coords[offset+1] = u[3] * x_tmp_0 + u[4] * x_tmp_1 + u[5] * x_tmp_2;
		coords[offset+2] = u[6] * x_tmp_0 + u[7] * x_tmp_1 + u[8] * x_tmp_2;
	}
}
__device__ floating_point_type calcRMS(
		floating_point_type* x,
		floating_point_type* y,
		unsigned int num_atoms){
	double sum_res = 0.0;

	for(unsigned int i=0; i<num_atoms*3; ++i){
		sum_res += (x[i] - y[i]) * (x[i] - y[i]);
	}

	return sqrt(sum_res/num_atoms);
}

///////////////////////////////////////////////////////////////
/// \remarks
/// Given a trajectory with N conformations stored in 'all_coordinates', it calculates the rmsd between the conformation
/// with 'i' = 'base_conformation_id' and the conformations in the range [j+1,M). In this case it tries to do each of the 
/// single rmsd measures in one thread.
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
__global__ void calcRMSDOfOneVsFollowing(
									 floating_point_type* first_conformation_coords,
									 const int base_conformation_id,
									 floating_point_type* all_coordinates,
									 const int number_of_conformations,
									 const int atoms_per_conformation,
									 const int coordinates_per_conformation,
									 floating_point_type* rmsd){
	
	int second_conformation_id = base_conformation_id + 1 + (blockDim.x*blockIdx.x + threadIdx.x);

	while(second_conformation_id < number_of_conformations){
		floating_point_type* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);

		rmsd[second_conformation_id-(base_conformation_id+1)] = calcRMSDOfTwoConformations(
				first_conformation_coords,
				second_conformation_coords,
				atoms_per_conformation);

		second_conformation_id += blockDim.x*gridDim.x;
	}
}



__global__ void calcRMSDOfOneVsFollowingWithRotation(
									 floating_point_type* first_conformation_coords,
									 const int base_conformation_id,
									 floating_point_type* all_coordinates,
									 const int number_of_conformations,
									 const int atoms_per_conformation,
									 const int coordinates_per_conformation,
									 floating_point_type* rmsd){

	int second_conformation_id = base_conformation_id + 1 + (blockDim.x*blockIdx.x + threadIdx.x);

	while(second_conformation_id < number_of_conformations){
		floating_point_type* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);
		floating_point_type rot_matrix[9];

		rmsd[second_conformation_id-(base_conformation_id+1)] = calcRMSDOfTwoConformations(
				first_conformation_coords,
				second_conformation_coords,
				atoms_per_conformation,
				rot_matrix);

		rotate3D(
				atoms_per_conformation,
				second_conformation_coords,
				rot_matrix);

		second_conformation_id += blockDim.x*gridDim.x;
	}
}

__global__ void calcRMSDOfOneVsFollowingFitDiffersCalc(
		floating_point_type* fitReference,
		floating_point_type* calcReference,
		int reference_conformation_number,
		floating_point_type* rmsd,
		int numberOfConformations,
		int coordinatesPerConformation,
		int atomsPerConformation,
		floating_point_type *allCoordinates,
		int coordinatesPerRMSDConformation,
		int atomsPerRMSDConformation,
		floating_point_type *allRMSDCoordinates){

	int second_conformation_id = reference_conformation_number + 1 + (blockDim.x*blockIdx.x + threadIdx.x);

	while(second_conformation_id < numberOfConformations){
		floating_point_type* second_conformation_coords = &(allCoordinates[second_conformation_id*coordinatesPerConformation]);
		floating_point_type rot_matrix[9];

		calcRMSDOfTwoConformations(
						fitReference,
						second_conformation_coords,
						atomsPerConformation,
						rot_matrix);

		floating_point_type* second_calc_conformation_coords = &(allRMSDCoordinates[second_conformation_id*coordinatesPerRMSDConformation]);

		rotate3D(
				atomsPerConformation,
				second_conformation_coords,
				rot_matrix);

		rotate3D(
				atomsPerRMSDConformation,
				second_calc_conformation_coords,
				rot_matrix);

		rmsd[second_conformation_id-(reference_conformation_number+1)] = calcRMS(
				calcReference,
				second_calc_conformation_coords,
				atomsPerRMSDConformation);

		second_conformation_id += blockDim.x*gridDim.x;
	}
}

