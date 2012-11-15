#include "statistics.h"

StatisticsCalculator::Statistics(double* data, int len) {
	data_size = len;
	data_array = data;
	mean=0;
    med=0;
    mode=0;
    var=0;
    skew=0;
    kurt=0;
}

StatisticsCalculator::~Statistics(){}

void StatisticsCalculator::calculateStatistics() {
	//	calculate (population-level) mean
	for (int i=0; i< count; ++i) {
        mean += array[i];
    }
    mean/=double(count);

    //	calculate (population-level) variance
	for (int i=0; i<count; ++i) {
		var+=pow((array[i]-mean), 2.0);
	}
	var/=count;

	//	calculate (population-level) skewness
	for (int i=0; i<count; ++i) {
		skew+=(array[i]-mean)*(array[i]-mean)*(array[i]-mean);
	}
	skew*=sqrt(double(count))/pow(var*double(count), 1.5);

	//	calculate (population-level) kurtosis
	for (int i=0; i<count; ++i) {
		kurt+=pow((array[i]-mean), 4.0);
	}
	kurt/=double(count)*pow(var, 2.0);
	kurt-=3.0;
}
