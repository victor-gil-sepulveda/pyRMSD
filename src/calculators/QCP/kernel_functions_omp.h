#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#include "kernel_functions_serial.h"

class ThRMSDSerialOmpKernel: public ThRMSDSerialKernel{
	public:
		ThRMSDSerialOmpKernel(){}
		virtual ~ThRMSDSerialOmpKernel(){}

		void calcRMSDOfOneVsFollowing(double* all_coordinates, int base_conformation_id,
				int other_conformations_starting_id, int number_of_conformations, int number_of_atoms,
				double* rmsd, int omp_threads);

};
#endif
