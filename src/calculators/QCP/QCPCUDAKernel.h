#ifndef __QCP_KERNEL_CUDA_H_
#define __QCP_KERNEL_CUDA_H_

#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __host__
#define __shared__
#define CUDA_KERNEL_DIM(...)

#else
#define CUDA_KERNEL_DIM(...)  <<< __VA_ARGS__ >>>

#endif


#include "../KernelFunctions.h"
#include <cstdlib>
#include <iostream>

#ifdef CUDA_PRECISION_SINGLE
	#define floating_point_type float
#else
	#define floating_point_type double
#endif

class QCPCUDAKernel: public KernelFunctions{

	public:
		QCPCUDAKernel(RMSDCalculationData* data,
						int threads_per_block,
						int blocks_per_grid);

		virtual ~QCPCUDAKernel();

		virtual void setCalculationCoords(RMSDCalculationData* data);


		void updateDeviceCoordinates(
				double * coordinates,
				floating_point_type* device_coordinates,
				float* tmp_coordinates,
				int coordinates_per_conformation,
				int number_of_conformations	);

		void updateHostCoordinates(
				double * coordinates,
				floating_point_type * device_coordinates,
				float * tmp_coordinates,
				int number_of_conformations,
				int coordinates_per_conformation
				);

		void updateHostRMSDs(
				int numberOfConformations,
				int reference_conformation_number,
				double* rmsd
				);

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

		virtual void matrixInit(RMSDCalculationData* data);

		virtual void matrixOneVsFollowingFitEqualCalc(
													double* reference,
													int reference_conformation_number,
													double* rmsd,
													RMSDCalculationData* data);

		virtual void matrixOneVsFollowingFitDiffersCalc(
														double* fitReference,
														double* calcReference,
														int reference_conformation_number,
														double* rmsd,
														RMSDCalculationData* data);
		int threads_per_block;
		int blocks_per_grid;

	private:
		QCPCUDAKernel(){}

	protected:
		float* tmpHostCoords;
		float* tmpHostRMSDs;

		floating_point_type* deviceCoords;
		floating_point_type* deviceRMSDs;
		floating_point_type* deviceReference;

		float* tmpCalcHostCoords;

		floating_point_type* deviceCalcCoords;
		floating_point_type* deviceCalcReference;

};


#endif
