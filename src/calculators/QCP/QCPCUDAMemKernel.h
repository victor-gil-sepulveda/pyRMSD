/*
 * QCPCUDAMemKernel.h
 *
 *  Created on: Apr 13, 2013
 *      Author: victor
 */

#ifndef QCPCUDAMEMKERNEL_H_
#define QCPCUDAMEMKERNEL_H_

#include "QCPCUDAKernel.h"

class QCPCUDAMemKernel: public QCPCUDAKernel {
	public:
		QCPCUDAMemKernel(
				RMSDCalculationData* data,
				int threads_per_block,
				int blocks_per_grid);

		virtual ~QCPCUDAMemKernel();

		void matrixInit(RMSDCalculationData* data);

		void matrixEnd(int , std::vector<double>& );

		void matrixOneVsFollowingFitEqualCalc(double* reference,
													int reference_conformation_number,
													double* rmsd,
													RMSDCalculationData* data);

		void matrixOneVsFollowingFitDiffersCalc(double* fitReference,
													   double* calcReference,
													   int reference_conformation_number,
													   double* rmsd,
													   RMSDCalculationData* data);
		floating_point_type* allDeviceRMSDs;
};

#endif /* QCPCUDAMEMKERNEL_H_ */
