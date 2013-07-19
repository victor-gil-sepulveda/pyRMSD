#ifndef _KERNEL_SERIAL_OMP_H_
#define _KERNEL_SERIAL_OMP_H_

#include "QCPSerialKernel.h"

class QCPOmpKernel: public QCPSerialKernel{
	public:
		QCPOmpKernel(int number_of_threads);
		virtual ~QCPOmpKernel(){}

		virtual void oneVsFollowingFitEqualCalcCoords(
					double* reference,
					int reference_conformation_number,
					double* rmsd,
					RMSDCalculationData* data);

		virtual void oneVsFollowingFitDiffersCalcCoords(
					double* fitReference,
					double* calcReference,
					int reference_conformation_number,
					double* rmsd,
					RMSDCalculationData* data);

		int number_of_threads;

	private:
		QCPOmpKernel(){number_of_threads =0;};
};
#endif
