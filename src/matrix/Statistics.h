#ifndef _STATISTICS_H_
#define _STATISTICS_H_

class StatisticsCalculator {
	public:
		StatisticsCalculator(float* data, int len);
		~StatisticsCalculator();
		void calculateStatistics();
		int data_size;
		float* data_array;
		double mean;
		double variance;
		double skewness;
		double kurtosis;
		float max;
		float min;
};

#endif
