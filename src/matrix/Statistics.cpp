/*
 *
 * Based on the work of  jjhaag AT dreamincode.net (2007)
 *
 *
 */

#include "Statistics.h"
#include <cmath>
#include <algorithm>

StatisticsCalculator::StatisticsCalculator(float* data, int len) {
	data_size = len;
	data_array = data;
	mean = 0;
	variance = 0;
    skewness = 0;
    kurtosis = 0;
}

StatisticsCalculator::~StatisticsCalculator(){}

void StatisticsCalculator::calculateStatistics() {
	//	calculate (population-level) mean
	for (int i=0; i< data_size; ++i) {
        mean += data_array[i];
    }
    mean/=double(data_size);

    //	calculate (population-level) variance
	for (int i=0; i<data_size; ++i) {
		variance+=pow((data_array[i]-mean), 2.0);
	}
	variance/=data_size;

	//	calculate (population-level) skewness
	for (int i=0; i<data_size; ++i) {
		skewness+=(data_array[i]-mean)*(data_array[i]-mean)*(data_array[i]-mean);
	}
	skewness*=sqrt(double(data_size))/pow(variance*double(data_size), 1.5);

	//	calculate (population-level) kurtosis
	for (int i=0; i<data_size; ++i) {
		kurtosis+=pow((data_array[i]-mean), 4.0);
	}
	kurtosis/=double(data_size)*pow(variance, 2.0);
	kurtosis-=3.0;

	max = *(std::max_element(data_array, data_array+data_size));
	min = *(std::min_element(data_array, data_array+data_size));
}
