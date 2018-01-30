#ifndef MATH_BASIC_VECTOR_MATH_H_
#define MATH_BASIC_VECTOR_MATH_H_

#include <limits>

#include "../containers/primitive_container.h"

namespace tachyon{
namespace math{

// Primitive container math
// min, max, mean, median
template <class T>
inline T min(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	T min_observed = container[0];
	for(U32 i = 1; i < container.size(); ++i)
		if(container[i] < min_observed) min_observed = container[i];

	return(min_observed);
}

template <class T>
inline T max(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	T max_observed = container[0];
	for(U32 i = 1; i < container.size(); ++i)
		if(container[i] > max_observed) max_observed = container[i];

	return(max_observed);
}

template <class T>
inline double mean(const containers::PrimitiveContainer<T>& container){
	if(container.size() == 0)
		return(container[0]);

	double total_observed = container[0];
	for(U32 i = 1; i < container.size(); ++i)
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

}
}

#endif /* MATH_BASIC_VECTOR_MATH_H_ */
