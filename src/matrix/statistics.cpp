#include "statistics.h"

//	constructor, with default values of an empty array and 0 elements
Statistics::Statistics(double* inputArray=NULL, long n=0) {
    if (n) {
        array=new double[n];
        count=n;
        for (int i=0; i<n; ++i) {
            array[i]=inputArray[i];
        }
    } else {
        array=NULL;
        count=0;
    }
    mean=0;
    med=0;
    mode=0;
    var=0;
    skew=0;
    kurt=0;
    sorted=false;
    meanCalculated=false;
    varCalculated=false;
}
;

//	destructor
Statistics::~Statistics() {
    this->count=0;
    delete[] this->array;
}
;

//	set data elements, with default values of an empty array and 0 elements
//	resets all variables to their default values, so any previous calculations
//	will be lost if new data is assigned to the object
void Statistics::setData(double* inputArray, long n) {
    delete[] array;
    if (n) {
        this->array=new double[n];
        this->count=n;
        for (int i=0; i<n; ++i) {
            this->array[i]=inputArray[i];
        }
    } else {
        array=NULL;
        this->count=0;
    }
    mean=0;
    med=0;
    mode=0;
    var=0;
    skew=0;
    kurt=0;
    sorted=false;
    meanCalculated=false;
    varCalculated=false;
}

//	calculate mean and store internally
void Statistics::setMean() {
    for (int i=0; i<count; ++i) {
        mean += array[i];
    }
    mean/=double(count);
    meanCalculated=true;
}
;

//	calculate median and store internally
void Statistics::setMedian() {
    if (!sorted) {
        std::sort(array, array+count);
        sorted=true;
    }

    if (count%2==1) {
        //	if odd number of elements, median is in the middle of the sorted array
        med=array[count/2];
    } else {
        //	if even number of elements, median is avg of two elements straddling
        //	the middle of the array
        med=(array[count/2] + array[(count/2)-1])/2;
    }
}
;

//	calculate mode and store internally
void Statistics::setMode() {
    if (!sorted) {
        std::sort(array, array+count);
        sorted=true;
    }
    std::deque<long> histCounts(1, 1);
    std::deque<double> histValues(1, array[0]);

    for (int i=1; i<count; ++i) {
        if (array[i] != histValues.back()) {
            histValues.push_back(array[i]);
            histCounts.push_back(1);
        } else {
            histCounts.back()++;
        }
    }

    long maxIndex=0;
    long maxCount=histCounts.at(0);
    for (int i=1; i<(long)histCounts.size(); ++i) {
        if (maxCount < histCounts.at(i)) {
            maxCount=histCounts.at(i);
            maxIndex=i;
        }
    }
    mode=histValues.at(maxIndex);
}
;

//	calculate (population-level) variance and store internally
void Statistics::setVariance() {
    if (!meanCalculated) {
        setMean();
    }
    for (int i=0; i<count; ++i) {
        var+=pow((array[i]-mean), 2.0);
    }
    var/=count;
}

//	calculate (population-level) skewness and store internally
void Statistics::setSkewness() {
    if (!meanCalculated) {
        setMean();
    }
    if (!varCalculated) {
        setVariance();
    }
    for (int i=0; i<count; ++i) {
        skew+=(array[i]-mean)*(array[i]-mean)*(array[i]-mean);
    }
    skew*=sqrt(double(count))/pow(var*double(count), 1.5);
}

//	calculate (population-level) kurtosis and store internally
void Statistics::setKurtosis() {
    if (!meanCalculated) {
        setMean();
    }
    if (!varCalculated) {
        setVariance();
    }
    for (int i=0; i<count; ++i) {
        kurt+=pow((array[i]-mean), 4.0);
    }
    kurt/=double(count)*pow(var, 2.0);
    kurt-=3.0;
}

//  calculate all statistics
void Statistics::calculateStatistics() {
    setMean();
    setMedian();
    setMode();
    setVariance();
    setSkewness();
    setKurtosis();
}
