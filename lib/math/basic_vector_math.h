#ifndef MATH_BASIC_VECTOR_MATH_H_
#define MATH_BASIC_VECTOR_MATH_H_

#include <limits>

#include "summary_statistics.h"
#include "containers/primitive_container.h"

namespace tachyon{
namespace math{

// Primitive container math
// min, max, mean, median
template <class T>
inline T min(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	T min_observed = container[0];
	for(uint32_t i = 1; i < container.size(); ++i)
		if(container[i] < min_observed) min_observed = container[i];

	return(min_observed);
}

template <class T>
inline T max(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	T max_observed = container[0];
	for(uint32_t i = 1; i < container.size(); ++i)
		if(container[i] > max_observed) max_observed = container[i];

	return(max_observed);
}

template <class T>
inline double mean(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	double total_observed = container[0];
	for(uint32_t i = 1; i < container.size(); ++i)
		total_observed += container[i];

	return(total_observed / container.size());
}

template <class T>
inline double average(const containers::PrimitiveContainer<T>& container){ return(mean(container)); }

template <class T>
inline T median(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	return(container[container.size()/2]);
}

// Todo: need to skip missing values and EOV values
template <class T>
inline SummaryStatistics summary_statistics(const containers::PrimitiveContainer<T>& container){
	math::SummaryStatistics ss;
	for(uint32_t i = 0; i < container.size(); ++i)
		ss += container[i];

	ss.calculate();
	return(ss);
}

template <class T>
inline SummaryStatistics summary_statistics(const containers::InfoContainer<T>& info_container){
	math::SummaryStatistics ss;
	for(uint32_t i = 0; i < info_container.size(); ++i){
		for(uint32_t j = 0; j < info_container[i].size(); ++j)
			ss += info_container[i][j];
	}

	ss.calculate();
	return(ss);
}

template <class T>
inline SummaryStatistics summary_statistics(const containers::FormatContainer<T>& format_container){
	math::SummaryStatistics ss;
	for(uint32_t i = 0; i < format_container.size(); ++i){
		for(uint32_t j = 0; j < format_container[i].size(); ++j){
			for(uint32_t k = 0; k < format_container[i][j].size(); ++k){
				ss += format_container[i][j][k];
			}
		}
	}

	ss.calculate();
	return(ss);
}

}
}

#endif /* MATH_BASIC_VECTOR_MATH_H_ */
