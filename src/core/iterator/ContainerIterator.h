#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

class ContainerIterator{
	typedef ContainerIterator self_type;
	typedef StreamContainer container_type;

	inline const bool isUniform(void) const{ return(this->container.header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container.header_stride.controller.uniform); }

public:
	U32 n_entries;
	U32 position;
	container_type& container;
};

}
}
}


#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
