/*
 * RMSDCalculatorTypes.h
 *
 *  Created on: 06/03/2013
 *      Author: victor
 */

#ifndef RMSDCALCULATORTYPES_H_
#define RMSDCALCULATORTYPES_H_
#include <string>

enum RMSDCalculatorType{
	KABSCH_SERIAL_CALCULATOR = -1,
	KABSCH_OMP_CALCULATOR,
	KABSCH_CUDA_CALCULATOR,

	QTRFIT_SERIAL_CALCULATOR,
	QTRFIT_OMP_CALCULATOR,
	QTRFIT_CUDA_CALCULATOR,

	QCP_SERIAL_CALCULATOR,
	QCP_OMP_CALCULATOR,
	QCP_CUDA_CALCULATOR
};

inline std::string calculatorTypeToString(RMSDCalculatorType type){
	switch (type) {
		case KABSCH_SERIAL_CALCULATOR:
			return "KABSCH_SERIAL_CALCULATOR";
			break;

		case KABSCH_OMP_CALCULATOR:
			return "KABSCH_OMP_CALCULATOR";
			break;

		case KABSCH_CUDA_CALCULATOR:
			return "KABSCH_CUDA_CALCULATOR";
			break;

		case QTRFIT_SERIAL_CALCULATOR:
			return "QTRFIT_SERIAL_CALCULATOR";
			break;

		case QTRFIT_OMP_CALCULATOR:
			return "QTRFIT_OMP_CALCULATOR";
			break;

		case QTRFIT_CUDA_CALCULATOR:
			return "QTRFIT_CUDA_CALCULATOR";
			break;

		case QCP_SERIAL_CALCULATOR:
			return "QCP_SERIAL_CALCULATOR";
			break;

		case QCP_OMP_CALCULATOR:
			return "QCP_OMP_CALCULATOR";
			break;

		case QCP_CUDA_CALCULATOR:
			return "QCP_CUDA_CALCULATOR";
			break;

		default:
			return "";
	}
}

#endif /* RMSDCALCULATORTYPES_H_ */
