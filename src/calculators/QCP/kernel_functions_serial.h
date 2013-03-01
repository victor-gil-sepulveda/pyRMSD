#ifndef _KERNEL_SERIAL_H_
#define _KERNEL_SERIAL_H_

#include "KernelFunctions.h"

class ThRMSDSerialKernel: public KernelFunctions{
	public:
		ThRMSDSerialKernel(){}
		virtual ~ThRMSDSerialKernel(){}

		void centerCoordsOfAllConformations(int number_of_conformations, int number_of_atoms, double* all_coordinates);

		void centerCoords( double* all_cooordinates, int number_of_atoms);

		double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms);

		double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		virtual void calcRMSDOfOneVsFollowing( double* all_coordinates, int base_conformation_id, int other_conformations_starting_id,
											int number_of_conformations, int number_of_atoms, double* rmsd);
};
#endif
