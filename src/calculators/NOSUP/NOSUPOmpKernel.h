/*
 * KABSCHSerialKernel.h
 *
 *  Created on: 14/10/2013
 *      Author: victor
 */

#ifndef NOSUPOMPKERNEL_H_
#define NOSUPOMPKERNEL_H_

#include "../KernelFunctions.h"

#include <iostream>
using namespace std;

class NoSuperpositionOmpKernel: public KernelFunctions {
	public:
		int number_of_threads;

		NoSuperpositionOmpKernel(int number_of_threads);
		virtual ~NoSuperpositionOmpKernel();

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
		void centerAllAtCOM(RMSDCalculationData*){cout<<"OMP only fit no center"<<endl;}
		void centerAllAtFittingCOM(RMSDCalculationData*){cout<<"OMP fit and calc no center"<<endl;}

};

#endif /* KABSCHSERIALKERNEL_H_ */
