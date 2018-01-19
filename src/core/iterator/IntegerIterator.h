#ifndef ITERATOR_INTEGERITERATOR_H_
#define ITERATOR_INTEGERITERATOR_H_

namespace Tachyon{
namespace Iterator{

class IntegerIterator{
private:
	typedef IntegerIterator self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::Container container_type;

	typedef const U64 (self_type::*currentFunctionType)(void) const;

public:
	IntegerIterator() :
		current_position(0),
		current_byte_position(0),
		l_buffer(0),
		current_type_width(4),
		current_type(Core::YON_TYPE_32B),
		buffer_data(nullptr),
		currentFunction(&self_type::__current<U32>)
	{}

	IntegerIterator(const char* const buffer, const U64 l_buffer) :
		current_position(0),
		current_byte_position(0),
		l_buffer(l_buffer),
		current_type_width(4),
		current_type(Core::YON_TYPE_32B),
		buffer_data(buffer),
		currentFunction(&self_type::__current<U32>)
	{}

	IntegerIterator(const buffer_type& buffer) :
		current_position(0),
		current_byte_position(0),
		l_buffer(buffer.pointer),
		current_type_width(4),
		current_type(Core::YON_TYPE_32B),
		buffer_data(buffer.data),
		currentFunction(&self_type::__current<U32>)
	{}

	IntegerIterator(const container_type& container) :
		current_position(0),
		current_byte_position(0),
		l_buffer(container.buffer_data_uncompressed.pointer),
		current_type_width(4),
		current_type(container.getDataPrimitiveType()),
		buffer_data(container.buffer_data_uncompressed.data),
		currentFunction(nullptr)
	{
		switch(container.getDataPrimitiveType()){
		case(Core::YON_TYPE_CHAR):
		case(Core::YON_TYPE_8B):  this->currentFunction = &self_type::__current<BYTE>; this->current_type_width = 1; break;
		case(Core::YON_TYPE_16B): this->currentFunction = &self_type::__current<U16>;  this->current_type_width = 2; break;
		case(Core::YON_TYPE_32B): this->currentFunction = &self_type::__current<U32>;  this->current_type_width = 4; break;
		case(Core::YON_TYPE_64B): this->currentFunction = &self_type::__current<U64>;  this->current_type_width = 8; break;
		default: this->currentFunction = &self_type::__current<U32>;  this->current_type_width = 4; break;
		}
	}

	~IntegerIterator(){}

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
		case(Core::YON_TYPE_8B):  this->currentFunction = &self_type::__current<BYTE>; this->current_type_width = 1; break;
		case(Core::YON_TYPE_16B): this->currentFunction = &self_type::__current<U16>;  this->current_type_width = 2; break;
		case(Core::YON_TYPE_32B): this->currentFunction = &self_type::__current<U32>;  this->current_type_width = 4; break;
		case(Core::YON_TYPE_64B): this->currentFunction = &self_type::__current<U64>;  this->current_type_width = 8; break;
		default: return(false);
		}
		return(true);
	}

	inline void increment(const U32& p){
		if(this->current_byte_position + p > this->l_buffer)
			return;

		++this->current_position;
		this->current_byte_position += p;
	}

	void operator++(void){
		if(this->current_byte_position + this->current_type_width > this->l_buffer)
			return;

		++this->current_position;
		this->current_byte_position += this->current_type_width;
	}

	inline const U64 size(void) const{ return(this->l_buffer); }

private:
	template <class T>
	inline const U64 __current(void) const{
		return(*reinterpret_cast<const T* const>(&this->buffer_data[this->current_byte_position]));
	}

private:
	U32 current_position;
	U64 current_byte_position;
	U64 l_buffer;
	BYTE current_type_width;
	Core::TACHYON_CORE_TYPE current_type;
	const char* const buffer_data;
	currentFunctionType currentFunction;
};

}
}



#endif /* ITERATOR_INTEGERITERATOR_H_ */
