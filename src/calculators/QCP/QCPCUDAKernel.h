#ifndef __QCP_KERNEL_CUDA_H_
#define __QCP_KERNEL_CUDA_H_

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
		QCPCUDAKernel(int threads_per_block, int blocks_per_grid){
			this->threads_per_block = threads_per_block;
			this->blocks_per_grid = blocks_per_grid;
			tmpHostCoords = deviceCoords = tmpHostRMSDs = deviceRMSDs = NULL;
		}

		virtual ~QCPCUDAKernel();

		virtual void init(
				double* coordinates,
				int atomsPerConformation,
				int coordinatesPerConformation,
				int numberOfConformations);

		virtual void changeCalculationCoords(double*);

		void updateDeviceCoordinates(
				int numberOfConformations,
				int coordinatesPerConformation,
				double * allCoordinates
				);

		void updateHostRMSDs(
				int numberOfConformations,
				int reference_conformation_number,
				double* rmsd
				);


		double calcRMSDOfTwoConformations( double* first_conformation_coords, double* second_conformation_coords,
				int number_of_atoms, double* rot_matrix = NULL);

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
			tmpHostCoords = deviceCoords = tmpHostRMSDs = deviceRMSDs = NULL;
		}

		floating_point_type* tmpHostCoords;
		floating_point_type* tmpHostRMSDs;

		floating_point_type* deviceCoords;
		floating_point_type* deviceRMSDs;
};
#endif
