#include "kernel_functions_cuda.cuh"

#define floating_point_type float

// Centers all conformations
__global__ void centerCoordsOfAllConformations(const int number_of_conformations, const int number_of_atoms, floating_point_type* all_coordinates){
	int current_conformation_id = blockDim.x*blockIdx.x + threadIdx.x;
	while(current_conformation_id < number_of_conformations){
		int coordinates_per_conformation = number_of_atoms * 3;
		floating_point_type* conformation_coordinates = &(all_coordinates[current_conformation_id*coordinates_per_conformation]);
		centerCoords(conformation_coordinates, number_of_atoms);
		current_conformation_id += blockDim.x*gridDim.x;
	}
}

// Centers ONE conformation given by coords[number_of_conformation]
// to be changed by coords[number_of_conformation*3*number_of_atoms_per_conformation]
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

// Calculates the rmsd of second over first conformations in 'innerProduct'
__device__ floating_point_type calcRMSDForTwoConformationsWithTheobaldMethod(floating_point_type *A, 
																			 const floating_point_type E0, 
																			 const int number_of_atoms){
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

    return sqrt(fabs(2.0f * (E0 - mxEigenV)/number_of_atoms));
}

__device__ floating_point_type calcRMSDOfTwoConformations(floating_point_type* first_conformation_coords, 
									 			floating_point_type* second_conformation_coords, 
									 			const int number_of_atoms){

	floating_point_type A[9];
	floating_point_type E0 = innerProduct(A, first_conformation_coords, second_conformation_coords, number_of_atoms);
	return calcRMSDForTwoConformationsWithTheobaldMethod(A, E0, number_of_atoms);
}

__global__ void calcRMSDOfOneVsOthers(floating_point_type* all_coordinates,
									 const int base_conformation_id,
									 const int other_conformations_starting_id,
									 const int number_of_conformations,
									 const int atoms_per_conformation,
									 const int coordinates_per_conformation,
									 floating_point_type* rmsd){
	
	// Usign shared mem for sructure coords doesn't help too much
	//extern __shared__ floating_point_type first_conformation_fast_coords[];
	/*if (threadIdx.x == 0){
		coords_per_conformation = coordinates_per_conformation;
		
		floating_point_type* first_conformation_coords = &(all_coordinates[base_conformation_id*coords_per_conformation]);
		for(int i = 0; i < coords_per_conformation; ++i){
			first_conformation_fast_coords[i] = first_conformation_coords[i];
		}
	}*/
	
	
	floating_point_type* first_conformation_coords = &(all_coordinates[base_conformation_id*coordinates_per_conformation]);
	int second_conformation_id = other_conformations_starting_id + (blockDim.x*blockIdx.x + threadIdx.x);
	while(second_conformation_id < number_of_conformations){
		floating_point_type* second_conformation_coords = &(all_coordinates[second_conformation_id*coordinates_per_conformation]);
		rmsd[second_conformation_id] = calcRMSDOfTwoConformations(first_conformation_coords, second_conformation_coords, atoms_per_conformation);
		second_conformation_id += blockDim.x*gridDim.x;
	}
}
