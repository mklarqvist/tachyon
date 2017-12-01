#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

class ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataInterface self_type;
	typedef StreamContainerHeader header_type;
	typedef IO::BasicBuffer buffer_type;

public:
	ContainerIteratorDataInterface();
	virtual ~ContainerIteratorDataInterface();

public:
	U32 position;        // iterator position
	U32 n_entries;       // size
	buffer_type& buffer; // buffer reference
};

class ContainerIteratorDataBase : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataBase self_type;

public:
	ContainerIteratorDataBase();
	virtual ~ContainerIteratorDataBase();

public:
	header_type& header; // header reference
};

template <class T = BYTE>
class ContainerIteratorData : public ContainerIteratorDataBase{
private:
	typedef ContainerIteratorData self_type;
	typedef const T* const_pointer;
	typedef const T* const const_pointer_final;
	typedef const T& const_reference;

public:
	ContainerIteratorData();
	~ContainerIteratorData();

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

// Todo: need stride iterator, recastable to correct return type
// Todo: need buffer iterator, recastable to correct return type

class ContainerIterator{
private:
	typedef ContainerIterator self_type;
	typedef ContainerIteratorDataBase data_iterator_type;

protected:
	typedef StreamContainer container_type;

public:
	inline const bool isUniform(void) const{ return(this->container.header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container.header_stride.controller.uniform); }
	inline const bool hasStrides(void) const{ return(this->container.header.controller.mixedStride); }

public:
	U32 n_entries;
	U32 position;
	container_type& container;
	data_iterator_type* data_iterator; // recast me as the correct base type
};

}
}
}

#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
