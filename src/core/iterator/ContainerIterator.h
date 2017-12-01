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
	ContainerIteratorDataInterface(buffer_type& buffer) : position(0), n_entries(0), buffer(buffer){}
	virtual ~ContainerIteratorDataInterface(){}

public:
	U32 position;        // iterator position
	U32 n_entries;       // size
	buffer_type& buffer; // buffer reference
};


template <class T>
class ContainerIteratorData : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorData self_type;
	typedef const T* const_pointer;
	typedef const T* const const_pointer_final;
	typedef const T& const_reference;

public:
	ContainerIteratorData(buffer_type& buffer) :
		ContainerIteratorDataInterface(buffer)
	{
		this->n_entries = buffer.pointer / sizeof(T);
		assert(buffer.pointer % sizeof(T) == 0);
	}
	~ContainerIteratorData(){}

	inline const_reference current(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position*sizeof(T)])); }
	inline const_pointer currentAt(void) const{ return(reinterpret_cast<const_pointer>(&this->buffer.data[this->position*sizeof(T)])); }
	inline const_reference first(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[0])); }
	inline const_reference last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[(this->n_entries - 1)*sizeof(T)]));
	}
	inline const_reference operator[](const U32& p) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[p*sizeof(T)])); }
	inline const_pointer   at(const U32& p) const{ return( reinterpret_cast<const_pointer>(&this->buffer.data[p*sizeof(T)])); }
	inline void operator++(){
		if(this->position + 1 == this->n_entries) return;
		++this->position;
	}

	inline void operator--(){
		if(this->position == 0) return;
		--this->position;
	}

	inline void operator+=(const U32& p){
		if(this->position + p >= this->n_entries){
			this->position = this->n_entries - 1;
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

template <class dataType, class strideType>
class ContainerIterator{
private:
	typedef ContainerIterator                 self_type;
	typedef ContainerIteratorData<dataType>   data_iterator_type;
	typedef ContainerIteratorData<strideType> stride_iterator_type;

protected:
	typedef StreamContainer container_type;

public:
	explicit ContainerIterator(void) :
		n_entries(0),
		position(0),
		fixed_stride_size(1),
		hasStrideIteratorSet(false),
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
		this->hasStrideIteratorSet = false;

		this->data_iterator = new data_iterator_type(this->container->buffer_data_uncompressed, this->container->header);

		// Construct this iterator if there is a mixed stride
		// otherwise we assume that the stride size is one (1)
		if(this->container->header.controller.mixedStride){
			this->hasStrideIteratorSet = true;
			// Stride iterator
			this->stride_iterator = new stride_iterator_type(this->container->buffer_strides_uncompressed, this->container->header_stride);
		}
	}

	inline const bool isUniform(void) const{ return(this->container->header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container->header_stride.controller.uniform); }
	inline const bool hasStrides(void) const{ return(this->container->header.controller.mixedStride); }

	inline const strideType& currentStride(void) const{
		// Assumingly uniform step size 1
		if(!this->hasStrideIteratorSet)
			return(this->fixed_stride_size);

		return(this->stride_iterator->current());
	}
	inline const dataType& currentData(void) const{ return(this->data_iterator->current()); }
	inline const dataType* currentDataStart(void) const{ return(this->data_iterator->currentAt()); }

	void operator++(void){
		++(*this->data_iterator);
	}

public:
	U32 n_entries; // number of entries
	U32 position;  // iterator position
	bool hasStrideIteratorSet;
	strideType fixed_stride_size;
	container_type* container;              // reference container
	data_iterator_type* data_iterator;      // recast me as the correct base type
	stride_iterator_type* stride_iterator;
};

}
}
}

#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
