#include "format_container_string.h"

namespace tachyon{
namespace containers{

FormatContainer<std::string>::FormatContainer() :
	__containers(nullptr)
{

}

FormatContainer<std::string>::FormatContainer(const data_container_type& data_container,
                                              const meta_container_type& meta_container,
                                                const std::vector<bool>& pattern_matches,
                                                              const uint64_t  n_samples) :
	__containers(nullptr)
{
	if(data_container.data_uncompressed.size() == 0)
		return;

	if(data_container.header.data_header.HasMixedStride()){
		this->__setupBalanced(data_container, meta_container, pattern_matches, n_samples);
	} else {
		this->__setupBalanced(data_container, meta_container, pattern_matches, n_samples, data_container.header.data_header.stride);
	}
}

FormatContainer<std::string>::FormatContainer(const data_container_type& container, const uint64_t n_samples) :
	__containers(nullptr)
{
	if(container.data_uncompressed.size() == 0)
		return;

	if(container.header.data_header.controller.mixedStride){
		this->__setup(container, n_samples);
	} else {
		this->__setup(container, n_samples, container.header.data_header.GetStride());
	}
}

FormatContainer<std::string>::~FormatContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__containers + i)->~PrimitiveGroupContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

void FormatContainer<std::string>::__setup(const data_container_type& container, const uint64_t& n_samples){
	if(container.strides_uncompressed.size() == 0)
		return;

	this->n_capacity = container.data_uncompressed.size() / n_samples;
	this->n_entries  = 0;

	if(this->n_capacity == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_capacity*sizeof(value_type)));
	stride_container_type strides(container);

	uint32_t current_offset = 0;
	uint32_t current_position = 0;
	while(true){
		//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << (this->*func)(container.buffer_strides_uncompressed, i) << std::endl;
		new( &this->__containers[current_position] ) value_type( container, current_offset, n_samples, strides[current_position] );
		current_offset += strides[current_position] * n_samples;
		++this->n_entries;
		if(current_offset == container.data_uncompressed.size()) break;
		assert(current_offset <= container.data_uncompressed.size());
		++current_position;
	}
	assert(current_offset == container.data_uncompressed.size());
}

void FormatContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint64_t& n_samples){
		this->n_entries = meta_container.size();
		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
		stride_container_type strides(data_container);

		uint32_t current_offset = 0;
		uint32_t strides_offset = 0;
		for(uint32_t i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].GetFormatPatternId()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, strides[strides_offset] );
				current_offset += strides[strides_offset] * n_samples;
				++strides_offset;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.data_uncompressed.size());
}

void FormatContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint64_t& n_samples, const uint32_t stride_size){
	this->n_entries = meta_container.size();
	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	uint32_t current_offset = 0;
	// Case 1: if data is uniform
	if(data_container.header.data_header.IsUniform()){
		for(uint32_t i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].GetFormatPatternId()]){
				new( &this->__containers[i] ) value_type( data_container, 0, n_samples, stride_size );
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}

		current_offset += stride_size * n_samples;
	}
	// Case 2: if data is not uniform
	else {
		for(uint32_t i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].GetFormatPatternId()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, stride_size );
				current_offset += stride_size * n_samples;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
	}
	assert(current_offset == data_container.data_uncompressed.size());
}

void FormatContainer<std::string>::__setup(const data_container_type& container, const uint64_t& n_samples, const uint32_t stride_size){
	this->n_entries = container.data_uncompressed.size() / n_samples / stride_size;

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries * sizeof(value_type)));

	uint32_t current_offset = 0;
	// Case 1: data is uniform -> give all samples the same value
	if(container.header.data_header.IsUniform()){
		for(uint32_t i = 0; i < this->n_entries; ++i)
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );

	}
	// Case 2: data is not uniform -> interpret data
	else {
		for(uint32_t i = 0; i < this->n_entries; ++i){
			//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << "fixed: " << stride_size << std::endl;
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );
			current_offset += stride_size * n_samples;
		}
	}
	assert(current_offset == container.data_uncompressed.size());
}

}
}
