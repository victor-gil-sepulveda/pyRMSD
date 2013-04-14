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
			double* coordinates,
			int atomsPerConformation,
			int coordinatesPerConformation,
			int numberOfConformations,
			int threads_per_block,
			int blocks_per_grid);

	virtual ~QCPCUDAMemKernel();

	void matrixInit(
					double* allFittingCoordinates,
					int coordinatesPerFittingConformation,
					double* allCalculationCoordinates,
					int coordinatesPerCalculationConformation,
					int numberOfConformations);

	void matrixEnd(double* rmsds_tmp, int rmsds_tmp_len, std::vector<double>& rmsds);

	void matrixOneVsFollowingFitEqualCalcWithoutConfRotation(
												double* reference,
												int reference_conformation_number,
												double* rmsd,
												int numberOfConformations,
												int coordinatesPerConformation,
												int atomsPerConformation,
												double *allCoordinates);

	void matrixOneVsFollowingFitDiffersCalcWithoutConfRotation(
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
	floating_point_type* allDeviceRMSDs;
};

#endif /* QCPCUDAMEMKERNEL_H_ */
