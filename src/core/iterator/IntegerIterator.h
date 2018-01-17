#ifndef ITERATOR_INTEGERITERATOR_H_
#define ITERATOR_INTEGERITERATOR_H_

namespace Tachyon{
namespace Iterator{

class IntegerIterator{
private:
	typedef IntegerIterator self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::StreamContainer container_type;

	typedef const U64 (self_type::*currentFunctionType)(void) const;

public:
	IntegerIterator();
	IntegerIterator(const char* const buffer, const U64 l_buffer);
	IntegerIterator(const buffer_type& buffer);
	IntegerIterator(const container_type& container);

	/**<
	 * Returns the current value as a U64 return type
	 * given the current byte stream.
	 * @return
	 */
	inline U64 current(void) const{ return((this->*currentFunction)()); }

	/**<
	 *
	 * @param primitive_type
	 * @return
	 */
	inline bool setType(const Core::TACHYON_CORE_TYPE primitive_type){
		this->current_type = primitive_type;
		switch(primitive_type){
		case(Core::YON_TYPE_CHAR):
		case(Core::YON_TYPE_8B):  this->currentFunction = &self_type::__current<BYTE>; break;
		case(Core::YON_TYPE_16B): this->currentFunction = &self_type::__current<U16>;  break;
		case(Core::YON_TYPE_32B): this->currentFunction = &self_type::__current<U32>;  break;
		case(Core::YON_TYPE_64B): this->currentFunction = &self_type::__current<U64>;  break;
		default: return(false);
		}
		return(true);
	}

	inline void increment(const U32& p){
		++this->current_position;
		this->current_byte_position += p;
	}

	void operator++(void);
	void operator--(void);

private:
	template <class T>
	inline const U64 __current(void) const{
		return(*reinterpret_cast<const T* const>(&this->buffer_data[this->current_byte_position]));
	}

private:
	U32 current_position;
	U64 current_byte_position;
	U64 l_buffer;
	Core::TACHYON_CORE_TYPE current_type;
	const char* const buffer_data;
	currentFunctionType currentFunction;
};

}
}



#endif /* ITERATOR_INTEGERITERATOR_H_ */
