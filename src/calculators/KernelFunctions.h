/*
 * KernelFunctions.h
 *
 *  Created on: 04/03/2013
 *      Author: victor
 */

#ifndef KERNELFUNCTIONS_H_
#define KERNELFUNCTIONS_H_
#include <vector>

class RMSDCalculationData;


class KernelFunctions{
	public:
		KernelFunctions(){}
		virtual ~KernelFunctions(){}

		virtual void setCalculationCoords(RMSDCalculationData* data){}

		virtual void oneVsFollowingFitEqualCalcCoords(
				double* reference,
				int reference_conformation_number,
				double* rmsd,
				RMSDCalculationData* data) = 0;

		virtual void oneVsFollowingFitDiffersCalcCoords(
				double* fitReference,
				double* calcReference,
				int reference_conformation_number,
				double* rmsd,
				RMSDCalculationData* data) = 0;

		virtual void handleSymmetriesWithCalcCoords(
				int reference_conformation_number,
				double* rmsd,
				RMSDCalculationData* data);

		virtual void matrixInit(RMSDCalculationData* data){}

		virtual void matrixEnd(int , std::vector<double>& ){}

		virtual void matrixOneVsFollowingFitEqualCalc(
													double* reference,
													int reference_conformation_number,
													double* rmsd,
													RMSDCalculationData* data){

			oneVsFollowingFitEqualCalcCoords(
									reference,
									reference_conformation_number,
									rmsd,
									data);
		}

		virtual void matrixOneVsFollowingFitDiffersCalc(
														double* fitReference,
														double* calcReference,
														int reference_conformation_number,
														double* rmsd,
														RMSDCalculationData* data){
			oneVsFollowingFitDiffersCalcCoords(
												fitReference,
												calcReference,
												reference_conformation_number,
												rmsd,
												data);
		}
		
		virtual void centerAllAtCOM(RMSDCalculationData*);
		virtual void centerAllAtFittingCOM(RMSDCalculationData*);
};


#endif /* KERNELFUNCTIONS_H_ */
