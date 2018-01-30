#ifndef MATH_SQUARE_MATRIX_H_
#define MATH_SQUARE_MATRIX_H_

#include <stdio.h> // size_t
#include "../support/type_definitions.h"

namespace tachyon{
namespace math{

template <class T>
class SquareMatrix{
private:
	typedef SquareMatrix self_type;

public:

	SquareMatrix(const U32 width) : __width(width), __data(new T*[__width])
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

}
}

#endif /* MATH_SQUARE_MATRIX_H_ */
