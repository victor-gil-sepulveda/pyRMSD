/*
 * KABSCHOmpKernel.h
 *
 *  Created on: 07/03/2013
 *      Author: victor
 */

#ifndef KABSCHOMPKERNEL_H_
#define KABSCHOMPKERNEL_H_

#include "../KABSCH/KABSCHSerialKernel.h"

class KABSCHOmpKernel: public KABSCHSerialKernel {
public:
	KABSCHOmpKernel(int number_of_threads);
	virtual ~KABSCHOmpKernel();


	virtual void init(
					double* coordinates,
					int atomsPerConformation,
					int coordinatesPerConformation,
					int numberOfConformations){}

	virtual void changeCalculationCoords(double*){}

	virtual void oneVsFollowingFitEqualCalcCoords(
			double* reference,
			int reference_conformation_number,
			double* rmsd,
			int numberOfConformations,
			int coordinatesPerConformation,
			int atomsPerConformation,
			double *allCoordinates);

	virtual void oneVsFollowingFitDiffersCalcCoords(
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
};

#endif /* KABSCHOMPKERNEL_H_ */
