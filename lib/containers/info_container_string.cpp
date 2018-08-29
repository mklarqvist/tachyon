#include "info_container_string.h"

namespace tachyon{
namespace containers{

InfoContainer<std::string>::InfoContainer() :
	containers_(nullptr)
{

}

InfoContainer<std::string>::InfoContainer(const data_container_type& container) :
	containers_(nullptr)
{
	if(container.header.data_header.HasMixedStride())
		this->Setup(container);
	else
		this->Setup(container, container.header.data_header.stride);
}

InfoContainer<std::string>::InfoContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches) :
	containers_(nullptr)
{
	if(data_container.header.data_header.HasMixedStride())
		this->SetupBalanced(data_container, meta_container, pattern_matches);
	else
		this->SetupBalanced(data_container, meta_container, pattern_matches, data_container.header.data_header.stride);
}

InfoContainer<std::string>::~InfoContainer(void){
for(std::size_t i = 0; i < this->n_entries; ++i)
	((this->containers_ + i)->~PrimitiveContainer)();

::operator delete[](static_cast<void*>(this->containers_));
}

// For mixed strides
void InfoContainer<std::string>::Setup(const data_container_type& container){
	if(container.strides_uncompressed.size() == 0)
		return;

	stride_container_type strides(container);
	this->n_entries = strides.size();

	if(this->size() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	uint32_t current_offset = 0;
	for(uint32_t i = 0; i < this->size(); ++i){
		new( &this->containers_[i] ) value_type(&container.data_uncompressed.data()[current_offset], strides[i]);
		current_offset += strides[i];
	}
	assert(current_offset == container.data_uncompressed.size());
}

void InfoContainer<std::string>::SetupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches){
	this->n_entries = meta_container.size();

	if(this->n_entries == 0)
		return;

	stride_container_type strides(data_container);
	this->containers_ = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	uint32_t current_offset = 0;
	uint32_t stride_offset = 0;

	for(uint32_t i = 0; i < this->size(); ++i){
		// Meta entry has no INFO
		if(meta_container[i].GetInfoPatternId() == -1){
			new( &this->containers_[i] ) value_type( );
		}
		// If pattern matches
		else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
			new( &this->containers_[i] ) value_type(&data_container.data_uncompressed.data()[current_offset], strides[stride_offset]);
			current_offset += strides[stride_offset];
			++stride_offset;
		}
		// Otherwise place an empty
		else {
			new( &this->containers_[i] ) value_type( );
		}
	}
	assert(current_offset == data_container.data_uncompressed.size());
}

// For fixed strides
void InfoContainer<std::string>::Setup(const data_container_type& container, const uint32_t stride_size){
	if(this->n_entries == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	uint32_t current_offset = 0;
	for(uint32_t i = 0; i < this->n_entries; ++i){
		//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
		new( &this->containers_[i] ) value_type(&container.data_uncompressed.data()[current_offset], stride_size);
		current_offset += stride_size;
	}
	assert(current_offset == container.data_uncompressed.size());
}

void InfoContainer<std::string>::SetupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint32_t stride_size){
	this->n_entries = meta_container.size();

	if(this->n_entries == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	if(data_container.header.data_header.IsUniform() == false){
		uint32_t current_offset = 0;
		for(uint32_t i = 0; i < this->n_entries; ++i){
			// If there are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->containers_[i] ) value_type( );
			} // If pattern matches
			else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->containers_[i] ) value_type(&data_container.data_uncompressed.data()[current_offset], stride_size);
				current_offset += stride_size;
			}
			// Otherwise place an empty
			else {
				new( &this->containers_[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.data_uncompressed.size());
	}
	// Data is uniform
	else {
		for(uint32_t i = 0; i < this->n_entries; ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->containers_[i] ) value_type(data_container.data_uncompressed.data(), stride_size);
			}
			// Otherwise place an empty
			else {
				new( &this->containers_[i] ) value_type( );
			}
		}
	}
}

}
}
