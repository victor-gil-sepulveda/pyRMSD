#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#include "kernel_functions_serial.h"

class ThRMSDSerialOmpKernel: public ThRMSDSerialKernel{
	public:
		ThRMSDSerialOmpKernel(){}
		virtual ~ThRMSDSerialOmpKernel(){}

		void calcRMSDOfOneVsFollowing( double* all_coordinates,
										   double* reference_conformation,
										   int reference_conformation_id,
										   int number_of_conformations,
										   int number_of_atoms,
										   double* rmsd);

};
#endif
