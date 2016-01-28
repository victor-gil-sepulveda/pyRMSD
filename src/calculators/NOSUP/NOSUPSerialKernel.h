/*
 * KABSCHSerialKernel.h
 *
 *  Created on: 14/10/2013
 *      Author: victor
 */

#ifndef NOSUPSERIALKERNEL_H_
#define NOSUPSERIALKERNEL_H_

#include "../KernelFunctions.h"

#include <iostream>
using namespace std;

class NoSuperpositionSerialKernel: public KernelFunctions {
	public:
		NoSuperpositionSerialKernel();
		virtual ~NoSuperpositionSerialKernel();

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
		
		// If not superposition is done, then centering must not be 
		// performed either
		void centerAllAtCOM(RMSDCalculationData*){cout<<"Serial only fit no center"<<endl;}
		void centerAllAtFittingCOM(RMSDCalculationData*){cout<<"Serial fit and calc no center"<<endl;}
};

#endif /* KABSCHSERIALKERNEL_H_ */
