#ifndef _KERNEL_SERIAL_H_
#define _KERNEL_SERIAL_H_

#define floating_point_type double
namespace ThRMSDSerialKernel{
	void centerCoordsOfAllConformations(int number_of_conformations, int number_of_atoms, floating_point_type* all_coordinates);
	void centerCoords( floating_point_type* all_cooordinates, int number_of_atoms);
	floating_point_type innerProduct(floating_point_type* A, floating_point_type* first_conformation_coords, floating_point_type* second_conformation_coords, int number_of_atoms);
	floating_point_type calcRMSDForTwoConformationsWithTheobaldMethod(floating_point_type *A, floating_point_type E0, int number_of_atoms);
	floating_point_type calcRMSDOfTwoConformations( floating_point_type* first_conformation_coords, floating_point_type* second_conformation_coords, int number_of_atoms);
	void calcRMSDOfOneVsOthers( floating_point_type* all_coordinates, int base_conformation_id, int other_conformations_starting_id,
											int number_of_conformations, int number_of_atoms, floating_point_type* rmsd);
}
#endif
