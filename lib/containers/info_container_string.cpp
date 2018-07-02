#include "info_container_string.h"

namespace tachyon{
namespace containers{

InfoContainer<std::string>::InfoContainer() :
	__containers(nullptr)
{

}

InfoContainer<std::string>::InfoContainer(const data_container_type& container) :
	__containers(nullptr)
{
if(container.header.data_header.hasMixedStride())
	this->__setup(container);
else
	this->__setup(container, container.header.data_header.stride);
}

InfoContainer<std::string>::InfoContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches) :
	__containers(nullptr)
{
if(data_container.header.data_header.hasMixedStride())
	this->__setupBalanced(data_container, meta_container, pattern_matches);
else
	this->__setupBalanced(data_container, meta_container, pattern_matches, data_container.header.data_header.stride);
}

InfoContainer<std::string>::~InfoContainer(void){
for(std::size_t i = 0; i < this->n_entries; ++i)
	((this->__containers + i)->~basic_string)();

::operator delete[](static_cast<void*>(this->__containers));
}

// For mixed strides
void InfoContainer<std::string>::__setup(const data_container_type& container){
	if(container.buffer_strides_uncompressed.size() == 0)
		return;

	stride_container_type strides(container);
	this->n_entries = strides.size();

	if(this->size() == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	U32 current_offset = 0;
	for(U32 i = 0; i < this->size(); ++i){
		new( &this->__containers[i] ) value_type(&container.buffer_data_uncompressed.data()[current_offset], strides[i]);
		current_offset += strides[i];
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

void InfoContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches){
	this->n_entries = meta_container.size();

	if(this->n_entries == 0)
		return;

	stride_container_type strides(data_container);
	this->__containers = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	U32 current_offset = 0;
	U32 stride_offset = 0;

	for(U32 i = 0; i < this->size(); ++i){
		// Meta entry has no INFO
		if(meta_container[i].getInfoPatternID() == -1){
			new( &this->__containers[i] ) value_type( );
		}
		// If pattern matches
		else if(pattern_matches[meta_container[i].getInfoPatternID()]){
			new( &this->__containers[i] ) value_type(&data_container.buffer_data_uncompressed.data()[current_offset], strides[stride_offset]);
			current_offset += strides[stride_offset];
			++stride_offset;
		}
		// Otherwise place an empty
		else {
			new( &this->__containers[i] ) value_type( );
		}
	}
	assert(current_offset == data_container.buffer_data_uncompressed.size());
}

// For fixed strides
void InfoContainer<std::string>::__setup(const data_container_type& container, const U32 stride_size){
	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	U32 current_offset = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
		new( &this->__containers[i] ) value_type(&container.buffer_data_uncompressed.data()[current_offset], stride_size);
		current_offset += stride_size;
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

void InfoContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U32 stride_size){
	this->n_entries = meta_container.size();

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	if(data_container.header.data_header.isUniform() == false){
		U32 current_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			// If there are no INFO fields
			if(meta_container[i].getInfoPatternID() == -1){
				new( &this->__containers[i] ) value_type( );
			} // If pattern matches
			else if(pattern_matches[meta_container[i].getInfoPatternID()]){
				new( &this->__containers[i] ) value_type(&data_container.buffer_data_uncompressed.data()[current_offset], stride_size);
				current_offset += stride_size;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.buffer_data_uncompressed.size());
	}
	// Data is uniform
	else {
		for(U32 i = 0; i < this->n_entries; ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].getInfoPatternID()]){
				new( &this->__containers[i] ) value_type(data_container.buffer_data_uncompressed.data(), stride_size);
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
	}
}

}
}
