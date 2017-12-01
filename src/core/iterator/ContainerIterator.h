#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

class ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataInterface self_type;

protected:
	typedef IO::BasicBuffer buffer_type;

public:
	ContainerIteratorDataInterface();
	virtual ~ContainerIteratorDataInterface();

public:
	U32 position;        // iterator position
	U32 n_entries;       // size
	buffer_type& buffer; // buffer reference
};

template <class T>
class ContainerIteratorBase : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorBase self_type;
	typedef const T* const_pointer;
	typedef const T* const const_pointer_final;
	typedef const T& const_reference;

public:
	ContainerIteratorBase();
	~ContainerIteratorBase();

	inline const_reference current(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position*sizeof(T)])); }
	inline const_reference first(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[0])); }
	inline const_reference last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[(this->n_entries - 1)*sizeof(T)]));
	}
	inline const_reference operator[](const U32& p) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[p*sizeof(T)])); }
	inline const_pointer   at(const U32& p) const{ return( reinterpret_cast<const_pointer>(&this->buffer.data[p*sizeof(T)])); }
	inline void operator++(){
		if(this->position == this->n_entries) return;
		++this->position;
	}

	inline void operator--(){
		if(this->position == 0) return;
		--this->position;
	}

	inline void operator+=(const U32& p){
		if(this->position + p > this->n_entries){
			this->position = this->n_entries;
			return;
		}
		this->position += p;
	}

	inline void operator-=(const U32& p){
		if(this->position - p < 0){
			this->position = 0;
			return;
		}
		this->position -= p;
	}
};

template <class T>
class ContainerIteratorData : public ContainerIteratorBase<T>{
private:
	typedef ContainerIteratorData self_type;
	typedef StreamContainerHeader header_type;

public:
	ContainerIteratorData(buffer_type& buffer, header_type& header);
	virtual ~ContainerIteratorData();

public:
	header_type& header; // header reference
};

template <class T>
class ContainerIteratorStride : public ContainerIteratorBase<T>{
private:
	typedef ContainerIteratorStride self_type;
	typedef StreamContainerHeaderStride header_type;

public:
	ContainerIteratorStride(buffer_type& buffer, header_type& header);
	virtual ~ContainerIteratorStride();

public:
	header_type& header; // header reference
};

class ContainerIterator{
private:
	typedef ContainerIterator self_type;
	typedef ContainerIteratorData<BYTE>   data_iterator_byte_type;
	typedef ContainerIteratorData<U16>    data_iterator_u16_type;
	typedef ContainerIteratorData<U32>    data_iterator_u32_type;
	typedef ContainerIteratorData<U64>    data_iterator_u64_type;
	typedef ContainerIteratorData<char>   data_iterator_char_type;
	typedef ContainerIteratorData<S16>    data_iterator_s16_type;
	typedef ContainerIteratorData<S32>    data_iterator_s32_type;
	typedef ContainerIteratorData<float>  data_iterator_float_type;
	typedef ContainerIteratorData<double> data_iterator_double_type;
	typedef ContainerIteratorStride<BYTE> stride_iterator_byte_type;
	typedef ContainerIteratorStride<U16>  stride_iterator_u16_type;
	typedef ContainerIteratorStride<U32>  stride_iterator_u32_type;
	typedef ContainerIteratorStride<U64>  stride_iterator_u64_type;
	typedef ContainerIteratorStride<char> stride_iterator_char_type;
	typedef ContainerIteratorStride<S16>  stride_iterator_s16_type;
	typedef ContainerIteratorStride<S32>  stride_iterator_s32_type;
	typedef ContainerIteratorStride<float>   stride_iterator_float_type;
	typedef ContainerIteratorStride<double>  stride_iterator_double_type;

protected:
	typedef StreamContainer container_type;

public:
	explicit ContainerIterator(void) :
		n_entries(0),
		position(0),
		container(nullptr),
		data_iterator(nullptr),
		stride_iterator(nullptr)
	{

	}
	~ContainerIterator(){
		// Do not delete container. It's allocated externally
		delete this->data_iterator;
		delete this->stride_iterator;
	}

	void operator()(container_type& c){
		// Store container reference
		this->container = &c;

		// Recycling is possible
		delete this->data_iterator;   this->data_iterator   = nullptr;
		delete this->stride_iterator; this->stride_iterator = nullptr;

		// Data iterator
		if(this->container->header.controller.signedness == 0){
			switch(this->container->header.controller.type){
			case(Core::TYPE_8B)  : this->data_iterator = new data_iterator_byte_type(this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_16B) : this->data_iterator = new data_iterator_u16_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_32B) : this->data_iterator = new data_iterator_u32_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_64B) : this->data_iterator = new data_iterator_u64_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_FLOAT)  : this->data_iterator = new data_iterator_float_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_DOUBLE) : this->data_iterator = new data_iterator_double_type(this->container->buffer_data_uncompressed, this->container->header); break;
			default: std::cerr << "illegal type!" << std::endl; exit(1); break;
			}
		} else {
			switch(this->container->header.controller.type){
			case(Core::TYPE_8B)  : this->data_iterator = new data_iterator_char_type(this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_16B) : this->data_iterator = new data_iterator_s16_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_32B) : this->data_iterator = new data_iterator_s32_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_FLOAT)  : this->data_iterator = new data_iterator_float_type (this->container->buffer_data_uncompressed, this->container->header); break;
			case(Core::TYPE_DOUBLE) : this->data_iterator = new data_iterator_double_type(this->container->buffer_data_uncompressed, this->container->header); break;
			default: std::cerr << "illegal type!" << std::endl; exit(1); break;
			}
		}

		// Construct this iterator if there is a mixed stride
		// otherwise we assume that the stride size is one (1)
		if(this->container->header.controller.mixedStride){
			// Stride iterator
			if(this->container->header_stride.controller.signedness == 0){
				switch(this->container->header_stride.controller.type){
				case(Core::TYPE_8B)  : this->stride_iterator = new stride_iterator_byte_type(this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_16B) : this->stride_iterator = new stride_iterator_u16_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_32B) : this->stride_iterator = new stride_iterator_u32_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_64B) : this->stride_iterator = new stride_iterator_u64_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_FLOAT)  : this->stride_iterator = new stride_iterator_float_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_DOUBLE) : this->stride_iterator = new stride_iterator_double_type(this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				default: std::cerr << "illegal type!" << std::endl; exit(1); break;
				}
			} else {
				switch(this->container->header_stride.controller.type){
				case(Core::TYPE_8B)  : this->stride_iterator = new stride_iterator_char_type(this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_16B) : this->stride_iterator = new stride_iterator_s16_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_32B) : this->stride_iterator = new stride_iterator_s32_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_FLOAT)  : this->stride_iterator = new stride_iterator_float_type (this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				case(Core::TYPE_DOUBLE) : this->stride_iterator = new stride_iterator_double_type(this->container->buffer_strides_uncompressed, this->container->header_stride); break;
				default: std::cerr << "illegal type!" << std::endl; exit(1); break;
				}
			}
		}
	}

	inline const bool isUniform(void) const{ return(this->container->header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container->header_stride.controller.uniform); }
	inline const bool hasStrides(void) const{ return(this->container->header.controller.mixedStride); }

public:
	U32 n_entries; // number of entries
	U32 position;  // iterator position
	container_type* container;                  // reference container
	data_iterator_byte_type* data_iterator;     // recast me as the correct base type
	stride_iterator_byte_type* stride_iterator; // recast me as the correct base type
};

}
}
}

#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
