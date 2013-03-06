#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#include "QCPSerialKernel.h"

class ThRMSDSerialOmpKernel: public ThRMSDSerialKernel{
	public:
		ThRMSDSerialOmpKernel(){}
		virtual ~ThRMSDSerialOmpKernel(){}

		virtual void oneVsFollowingFitEqualCalcWithoutConfRotation(
							double* reference,
							int reference_conformation_number,
							double* rmsd,
							int numberOfConformations,
							int coordinatesPerConformation,
							int atomsPerConformation,
							double *allCoordinates);

		virtual void oneVsFollowingFitEqualCalcWithConfRotation(
					double* reference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates);

		virtual void oneVsFollowingFitDiffersCalcWithoutConfRotation(
					double* fitReference,
					double* calcReference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates,
					int coordinatesPerRMSDConformation,
					int atomsPerRMSDConformation,
					double *allRMSDCoordinates);

		virtual void oneVsAllFitDiffersCalcWithConfRotation(
					double* fitReference,
					double* calcReference,
					int reference_conformation_number,
					double* rmsd,
					int numberOfConformations,
					int coordinatesPerConformation,
					int atomsPerConformation,
					double *allCoordinates,
					int coordinatesPerRMSDConformation,
					int atomsPerRMSDConformation,
					double *allRMSDCoordinates);

};
#endif
