#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

// Todo: need stride iterator, recastable to correct return type
// Todo: need buffer iterator, recastable to correct return type

class ContainerIterator{
private:
	typedef ContainerIterator self_type;
protected:
	typedef StreamContainer container_type;

public:
	inline const bool isUniform(void) const{ return(this->container.header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container.header_stride.controller.uniform); }
	inline const bool hasStrides(void) const{ return(this->container.header.controller.mixedStride); }

	void operator++();
	void operator--();

public:
	U32 n_entries;
	U32 position;
	container_type& container;
};

class ContainerIteratorBYTE : public ContainerIterator{
	typedef const BYTE& (self_type::*operator_function_type)(const U32& p);

	const BYTE& current(void) const;
	const BYTE& first(void) const;
	const BYTE& last(void) const;
	const BYTE& operator[](const U32& p) const;
	void operator++();
	void operator--();
};

}
}
}


#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
