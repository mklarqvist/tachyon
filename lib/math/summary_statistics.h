#ifndef SummaryStatistics_H_
#define SummaryStatistics_H_

#include <cmath>
#include <limits>

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

	template <class T>
	void operator+=(const T& value){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total += 1;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	template <class T>
	void add(const T& value, const double& weight = 1){
		this->total         += value;
		this->total_squared += value*value;
		this->n_total       += weight;
		if(value < this->min) this->min = value;
		if(value > this->max) this->max = value;
	}

	template <class T>
	void addNonzero(const T& value){
		if(value != 0){
			this->total         += value;
			this->total_squared += value*value;
			this->n_total       += 1;
			if(value < this->min) this->min = value;
			if(value > this->max) this->max = value;
		}
	}

	void reset(void){
		this->total         = 0;
		this->total_squared = 0;
		this->n_total       = 0;
		this->mean          = 0;
		this->standard_deviation = 0;
		this->min = std::numeric_limits<double>::max();
		this->max = std::numeric_limits<double>::min();
	}

	// Accessor functions
	inline double getTotal(void) const{ return(this->total); }
	inline double getTotalSquared(void) const{ return(this->total_squared); }
	inline double getCount(void) const{ return(this->n_total); }
	inline double getMean(void) const{ return(this->mean); }
	inline double getStandardDeviation(void) const{ return(this->standard_deviation); }
	inline double getMin(void) const{ return(this->min); }
	inline double getMax(void) const{ return(this->max); }

	friend std::ostream& operator<<(std::ostream& stream, self_type& self){
		self.calculate();
		stream << self.n_total << "\t" << self.mean << "\t" << self.standard_deviation << "\t" << self.min << "\t" << self.max;
		return(stream);
	}

public:
	double total;
	double total_squared;
	double n_total;
	double mean;
	double standard_deviation;
	double min;
	double max;
};

}
}



#endif /* SummaryStatistics_H_ */
