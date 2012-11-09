#ifndef _DEVICE_RMSD_2_H_
#define _DEVICE_RMSD_2_H_

#define floating_point_type float

__global__ void centerCoordsOfAllConformations(const int number_of_conformations, const int number_of_atoms, floating_point_type* all_coordinates);

__device__ void centerCoords( floating_point_type* all_cooordinates, const int number_of_atoms);

__device__ floating_point_type innerProduct(floating_point_type* A, 
											floating_point_type* first_conformation_coords, 
											floating_point_type* second_conformation_coords, 
											const int number_of_atoms);

__device__ floating_point_type calcRMSDForTwoConformationsWithTheobaldMethod(floating_point_type *A, 
																			const floating_point_type E0, 
																			const int number_of_atoms);

__device__ floating_point_type calcRMSDOfTwoConformations( floating_point_type* first_conformation_coords, 
														   floating_point_type* second_conformation_coords, 
														   const int number_of_atoms);

__global__ void calcRMSDOfOneVsOthers( 	floating_point_type* all_coordinates, const int base_conformation_id, 
										const int other_conformations_starting_id, 
										const int number_of_conformations, const int number_of_atoms, 
										const int atoms_per_conformation, floating_point_type* rmsd);
#endif