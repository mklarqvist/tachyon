#include "primitive_group_container_string.h"

namespace tachyon{
namespace containers{

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer() : __n_objects(0), __strings(nullptr){}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32& n_entries, const U32 strides_each) :
	__n_objects(n_entries), // limitation
	__strings(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
{
	U32 current_offset = offset;
	for(size_type i = 0; i < this->size(); ++i){
		new( &this->__strings[i] ) value_type( &container.buffer_data_uncompressed[current_offset], strides_each );
		current_offset += strides_each;
	}
}

PrimitiveGroupContainer<std::string>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->size(); ++i)
		((this->__strings + i)->~basic_string)();

	::operator delete[](static_cast<void*>(this->__strings));
}

}
}
