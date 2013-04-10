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
		QCPCUDAKernel(double* coordinates,
				int atomsPerConformation,
				int coordinatesPerConformation,
				int numberOfConformations,
				int threads_per_block,
				int blocks_per_grid);

		virtual ~QCPCUDAKernel();

		virtual void changeCalculationCoords(
				double* calcCoords,
				int number_of_atoms,
				int numberOfConformations);


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

		virtual void oneVsFollowingFitDiffersCalcWithConfRotation(
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

		double innerProduct(double* A, double* first_conformation_coords, double* second_conformation_coords, int number_of_atoms);

		double calcRMSDForTwoConformationsWithTheobaldMethod(double *A, double E0, int number_of_atoms, double* rot_matrix = NULL);

		int threads_per_block;
		int blocks_per_grid;

	private:

		QCPCUDAKernel(){
			threads_per_block = blocks_per_grid = 0;

			tmpHostCoords = deviceCoords =
			tmpHostRMSDs = deviceRMSDs =
			tmpHostReference = deviceReference =
			tmpCalcHostCoords = deviceCalcCoords =
			tmpCalcHostReference = deviceCalcReference = NULL;
		}

		floating_point_type* tmpHostCoords;
		floating_point_type* tmpHostRMSDs;
		floating_point_type* tmpHostReference;

		floating_point_type* deviceCoords;
		floating_point_type* deviceRMSDs;
		floating_point_type* deviceReference;

		floating_point_type* tmpCalcHostCoords;
		floating_point_type* tmpCalcHostReference;

		floating_point_type* deviceCalcCoords;
		floating_point_type* deviceCalcReference;
};
#endif
