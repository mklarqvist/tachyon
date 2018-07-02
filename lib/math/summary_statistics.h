#ifndef SummaryStatistics_H_
#define SummaryStatistics_H_

#include <cmath>
#include <limits>

#include "support/type_definitions.h"

namespace tachyon{
namespace math{

struct SummaryStatistics{
private:
	typedef SummaryStatistics self_type;

public:
	SummaryStatistics() :
		total(0),
		total_squared(0),
		n_total(0),
		mean(0),
		standard_deviation(0),
		min(std::numeric_limits<double>::max()),
		max(std::numeric_limits<double>::min())
	{

	}

	bool calculate(void){
		if(this->n_total == 0){
			this->mean = 0;
			this->standard_deviation = 0;
			return false;
		}

		this->mean = this->total / this->n_total;

		if(this->n_total > 1){
			this->standard_deviation = sqrt(this->total_squared/this->n_total - (this->total / this->n_total)*(this->total / this->n_total));
		} else this->standard_deviation = 0;

		return true;
	}

	inline double getSigma(void) const{ return(this->standard_deviation); }
	inline double getSigmaSquared(void) const{ return(this->standard_deviation*this->standard_deviation); }

	template <class T> void operator+=(const T& value){
		this->total         += value;
		this->total_squared += value*value;
		++this->n_total;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

public:
	double total;
	double total_squared;
	U64    n_total;
	double mean;
	double standard_deviation;
	double min;
	double max;
};

}
}



#endif /* SummaryStatistics_H_ */
