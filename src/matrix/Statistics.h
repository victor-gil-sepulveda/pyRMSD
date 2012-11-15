/******************************************************************************/
/* Statistics class
 *
 * jjhaag AT dreamincode.net
 *
 * 14 November, 2007
 *
 * Class for calculating basic statistics from an array of data elements.
 * Capable of calculating the three standard measures of central tendency (mean,
 * median, and mode), and three standard measures of the shape of a distribution
 * (variance, skewness, and kurtosis).
 *
 * The measures of shape are all population-level metrics, rather than sample
 * metrics.  This means that they will be biased for samples from a population.
 * Search for the appropriate terms on wikipedia or google to find how to change
 * these methods to return the sample variance, skewness, and kurtosis.
 *
 * And also note that the mode does not make sense for continuously-distributed
 * data; although the class accepts double-precision floating point numbers as
 * data, the mode will only make sense when the data fall into discrete classes.
 *
 *
 * The required data is an array of double-precision elements, and a count of
 * the number of elements in the array.
 *
 * For clarity, the class makes no checks as to whether the number of data
 * elements is sufficient to calculate each of the statistics, at least to
 * ensure that there are no infinite/indeterminate terms.  The measures of
 * central tendency only require a single element, though it's kind of pointless
 * in that case.  As they are presented here (i.e. as population-level
 * statistics), the measures of shape only require a single element.  However,
 * to convert to sample statistics, the variance would require a minimum of 2
 * elements, skewness would require at least 3 elements, and kurtosis would
 * require at least 4 elements.
 *
 * The class is sparse on comments within the code, but the methods and variable
 * names are hopefully fairly self-explanatory.  If this is not the case, leave
 * a comment and I will clarify (this is just in reference to the code itself;
 * if you are unclear on the statistical terminology, google it).
 *
 *
 * Dependencies:  requires <deque>, <algorithm>, and <cmath>.
 *
 * Example usage:
 *
 * double array[]={1.0,3.0,2.0,8.9,10.0,0.5};
 * Statistics stats(array, 8);
 * stats.calculateStatistics();
 */

#ifndef _STATISTICS_H_
#define _STATISTICS_H_


#include <deque>
#include <algorithm>
#include <cmath>

class StatisticsCalculator {
	public:
		Statistics(double* data, int len);
		~Statistics();
		void calculateStatistics();
		int data_size;
		double* data_array;
		double mean;
		double median;
		double mode;
		double variance;
		double skewness;
		double kurtosis;
}

#endif
