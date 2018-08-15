#include "primitive_group_container_string.h"

namespace tachyon{
namespace containers{

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer() : containers_(nullptr){}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const data_container_type& container,
	                                                          const U32& offset,
	                                                          const U32& n_entries,
	                                                          const U32 strides_each) :
	PrimitiveGroupContainerInterface(n_entries), // limitation
	containers_(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
{
	U32 current_offset = offset;
	for(size_type i = 0; i < this->size(); ++i){
		// check length
		size_type j = 0;
		for(; j < strides_each; ++j){
			// Find premature end-of-string marker
			if(container.buffer_data_uncompressed[current_offset + j] == '\0'){
				break;
			}
		}
		new( &this->containers_[i] ) value_type( &container.buffer_data_uncompressed[current_offset], j );
		current_offset += strides_each;
	}
}

PrimitiveGroupContainer<std::string>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->size(); ++i)
		((this->containers_ + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->containers_));
}

}
}
