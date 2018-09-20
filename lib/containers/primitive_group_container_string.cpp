#include "primitive_group_container_string.h"

namespace tachyon{
namespace containers{

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer() : containers_(nullptr){}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const self_type& other) :
	PrimitiveGroupContainerInterface(other),
	containers_(static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type))))
{
	for(int i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( other.at(i) );
}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const data_container_type& container,
	                                                          const uint32_t& offset,
	                                                          const uint32_t& n_entries,
	                                                          const uint32_t strides_each) :
	PrimitiveGroupContainerInterface(n_entries), // limitation
	containers_(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
{
	uint32_t current_offset = offset;
	for(size_type i = 0; i < this->size(); ++i){
		// check length
		size_type j = 0;
		for(; j < strides_each; ++j){
			// Find premature end-of-string marker and stop.
			if(container.data_uncompressed[current_offset + j] == '\0'){
				break;
			}
		}
		new( &this->containers_[i] ) value_type( &container.data_uncompressed[current_offset], j );
		current_offset += strides_each;
	}
}

PrimitiveGroupContainer<std::string>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->size(); ++i)
		((this->containers_ + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->containers_));
}

void PrimitiveGroupContainer<std::string>::resize(void){
	pointer temp       = this->containers_;
	this->n_capacity_ *= 2;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for(uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for(std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(temp));
}

void PrimitiveGroupContainer<std::string>::resize(const size_t new_size){
	// if new size < current capacity
	if(new_size < this->n_capacity_){
		// if new size < current number of entries
		if(new_size < this->n_objects_){
			this->n_objects_ = new_size;
			return;
		}
		return;
	}

	pointer temp       = this->containers_;
	this->n_capacity_  = new_size;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for(uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for(std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(temp));
}

}
}
