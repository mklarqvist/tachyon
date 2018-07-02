#ifndef MATH_SQUARE_MATRIX_H_
#define MATH_SQUARE_MATRIX_H_

#include <stdio.h> // size_t

#include "algorithm/permutation/permutation_manager.h"

namespace tachyon{
namespace math{

template <class T>
class SquareMatrix{
private:
	typedef SquareMatrix                  self_type;
	typedef algorithm::PermutationManager ppa_type;

public:
	SquareMatrix(const U32 width) :
		__width(width),
		__data(new T*[__width])
	{
		for(U32 i = 0; i < this->__width; ++i){
			this->__data[i] = new T[width];
			memset(this->__data[i], 0, sizeof(T)*this->__width);
		}
	}

	~SquareMatrix(){
		for(U32 i = 0; i < this->__width; ++i)
			delete [] this->__data[i];
		delete [] this->__data;
	}

	inline T& operator()(const U32& i, const U32& j){ return(this->__data[i][j]); }
	inline const T& operator()(const U32& i, const U32& j) const{ return(this->__data[i][j]); }

	// Basic math operators
	template <class Y>
	self_type& operator/=(const Y& value){
		if(value == 0)
			return(*this);

		for(U32 i = 0; i < this->__width; ++i){
			for(U32 j = 0; j < this->__width; ++j){
				this->__data[i][j] /= value;
			}
		}
		return(*this);
	}

	self_type& operator+=(const self_type& other);
	self_type& operator-=(const self_type& other);
	self_type& operator/=(const self_type& other);
	self_type& operator*=(const self_type& other);
	self_type& add(const self_type& other, const ppa_type& ppa_manager);
	self_type& addUpperTriagonal(const self_type& other, const ppa_type& ppa_manager);

	// Utility
	void clear(void){
		for(U32 i = 0; i < this->__width; ++i){
			for(U32 j = 0; j < this->__width; ++j)
				this->__data[i][j] = 0;
		}
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& matrix){
		for(U32 i = 0; i < matrix.__width; ++i){
			out << matrix.__data[i][0];
			for(U32 j = 1; j < matrix.__width; ++j){
				out << '\t' << matrix.__data[i][j];
			}
			out << '\n';
		}
		return(out);
	}

private:
	size_t __width;
	T**    __data;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
SquareMatrix<T>& SquareMatrix<T>::operator+=(const self_type& other){
	for(U32 i = 0; i < this->__width; ++i){
		for(U32 j = 0; j < this->__width; ++j){
			this->__data[i][j] += other.__data[i][j];
		}
	}
	return(*this);
}

template <class T>
SquareMatrix<T>& SquareMatrix<T>::operator-=(const self_type& other){
	for(U32 i = 0; i < this->__width; ++i){
		for(U32 j = 0; j < this->__width; ++j){
			this->__data[i][j] -= other.__data[i][j];
		}
	}
	return(*this);
}

template <class T>
SquareMatrix<T>& SquareMatrix<T>::operator/=(const self_type& other){
	for(U32 i = 0; i < this->__width; ++i){
		for(U32 j = 0; j < this->__width; ++j){
			this->__data[i][j] /= other.__data[i][j];
		}
	}
	return(*this);
}

template <class T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const self_type& other){
	for(U32 i = 0; i < this->__width; ++i){
		for(U32 j = 0; j < this->__width; ++j){
			this->__data[i][j] *= other.__data[i][j];
		}
	}
	return(*this);
}

template <class T>
SquareMatrix<T>& SquareMatrix<T>::add(const self_type& other, const ppa_type& ppa_manager){
	assert(ppa_manager.n_samples == this->__width);
	for(U32 i = 0; i < this->__width; ++i){
		for(U32 j = 0; j < this->__width; ++j){
			//std::cerr << "(" << i << "," << j << ") -> (" << ppa_manager[i] << "," << ppa_manager[j] << ")" << std::endl;
			this->__data[i][j] += other.__data[ppa_manager[i]][ppa_manager[j]];
		}
	}
	return(*this);
}

template <class T>
SquareMatrix<T>& SquareMatrix<T>::addUpperTriagonal(const self_type& other, const ppa_type& ppa_manager){
	assert(ppa_manager.n_samples == this->__width);
	for(U32 i = 0; i < this->__width; ++i){
		const U32& pA = ppa_manager[i];

		for(U32 j = i; j < this->__width; ++j){
			const U32& pB = ppa_manager[j];
			//std::cerr << i << "/" << j << '\t' << pA << '/' << pB << '\t' << (int)this->__data[i][j] << std::endl;
			if(pA >= pB) this->__data[pB][pA] += other.__data[i][j];
			else         this->__data[pA][pB] += other.__data[i][j];
		}
	}
	return(*this);
}

}
}

#endif /* MATH_SQUARE_MATRIX_H_ */
