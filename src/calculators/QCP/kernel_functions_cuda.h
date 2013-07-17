#ifndef _DEVICE_RMSD_2_H_
#define _DEVICE_RMSD_2_H_

#ifdef __CDT_PARSER__
	#define __global__
	#define __device__
	#define __host__
	#define __shared__
#endif


#ifdef CUDA_PRECISION_SINGLE
	#define floating_point_type float
#else
	#define floating_point_type double
#endif


__device__ floating_point_type innerProduct(
		floating_point_type* A,
		floating_point_type* first_conformation_coords,
		floating_point_type* second_conformation_coords,
		const int number_of_atoms);

__device__ floating_point_type calcRMSDForTwoConformationsWithTheobaldMethod(
		floating_point_type *A,
		const floating_point_type E0,
		const int number_of_atoms,
		floating_point_type* rot_matrix = NULL);

__device__ floating_point_type calcRMSDOfTwoConformations(
		floating_point_type* first_conformation_coords,
		floating_point_type* second_conformation_coords,
		const int number_of_atoms,
		floating_point_type* rot_matrix = NULL);

__device__ void rotate3D(
		unsigned int number_of_atoms,
		floating_point_type * const coords,
		floating_point_type* u);

__device__ floating_point_type calcRMS(
		floating_point_type* x,
		floating_point_type* y,
		unsigned int num_atoms);

__global__ void calcRMSDOfOneVsFollowing(
		 floating_point_type* first_conformation_coords,
		 const int base_conformation_id,
		 floating_point_type* all_coordinates,
		 const int number_of_conformations,
		 const int atoms_per_conformation,
		 const int coordinates_per_conformation,
		 floating_point_type* rmsd);

__global__ void calcRMSDOfOneVsFollowingWithRotation(
		 floating_point_type* first_conformation_coords,
		 const int base_conformation_id,
		 floating_point_type* all_coordinates,
		 const int number_of_conformations,
		 const int atoms_per_conformation,
		 const int coordinates_per_conformation,
		 floating_point_type* rmsd);

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
		floating_point_type *allRMSDCoordinates);
#endif
