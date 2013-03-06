#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#include "QCPSerialKernel.h"

class QCPOmpKernel: public QCPSerialKernel{
	public:
		QCPOmpKernel(int number_of_threads);
		virtual ~QCPOmpKernel(){}

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

		int number_of_threads;

	private:
		QCPOmpKernel(){number_of_threads =0;};
};
#endif
