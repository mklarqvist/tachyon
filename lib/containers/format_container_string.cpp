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
                                                              const U64  n_samples) :
	__containers(nullptr)
{
	if(data_container.buffer_data_uncompressed.size() == 0)
		return;

	if(data_container.header.data_header.hasMixedStride()){
		this->__setupBalanced(data_container, meta_container, pattern_matches, n_samples);
	} else {
		this->__setupBalanced(data_container, meta_container, pattern_matches, n_samples, data_container.header.data_header.stride);
	}
}

FormatContainer<std::string>::FormatContainer(const data_container_type& container, const U64 n_samples) :
	__containers(nullptr)
{
	if(container.buffer_data_uncompressed.size() == 0)
		return;

	if(container.header.data_header.controller.mixedStride){
		this->__setup(container, n_samples);
	} else {
		this->__setup(container, n_samples, container.header.data_header.getStride());
	}
}

FormatContainer<std::string>::~FormatContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__containers + i)->~PrimitiveGroupContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

void FormatContainer<std::string>::__setup(const data_container_type& container, const U64& n_samples){
	if(container.buffer_strides_uncompressed.size() == 0)
		return;

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
	stride_container_type strides(container);

	U32 current_offset = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << (this->*func)(container.buffer_strides_uncompressed, i) << std::endl;
		new( &this->__containers[i] ) value_type( container, current_offset, n_samples, strides[i] );
		current_offset += strides[i] * n_samples;
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

void FormatContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples){
		this->n_entries = meta_container.size();
		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
		stride_container_type strides(data_container);

		U32 current_offset = 0;
		U32 strides_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].getInfoPatternID() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].getFormatPatternID()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, strides[strides_offset] );
				current_offset += strides[strides_offset] * n_samples;
				++strides_offset;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.buffer_data_uncompressed.size());
}

void FormatContainer<std::string>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples, const U32 stride_size){
	this->n_entries = meta_container.size();
	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	U32 current_offset = 0;
	// Case 1: if data is uniform
	if(data_container.header.data_header.isUniform()){
		for(U32 i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].getInfoPatternID() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].getFormatPatternID()]){
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
		for(U32 i = 0; i < this->n_entries; ++i){
			// There are no INFO fields
			if(meta_container[i].getInfoPatternID() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].getFormatPatternID()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, stride_size );
				current_offset += stride_size * n_samples;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
	}
	assert(current_offset == data_container.buffer_data_uncompressed.size());
}

void FormatContainer<std::string>::__setup(const data_container_type& container, const U64& n_samples, const U32 stride_size){
	this->n_entries = container.buffer_data_uncompressed.size() / n_samples / stride_size;

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries * sizeof(value_type)));

	U32 current_offset = 0;
	// Case 1: data is uniform -> give all samples the same value
	if(container.header.data_header.isUniform()){
		for(U32 i = 0; i < this->n_entries; ++i)
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );

	}
	// Case 2: data is not uniform -> interpret data
	else {
		for(U32 i = 0; i < this->n_entries; ++i){
			//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << "fixed: " << stride_size << std::endl;
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );
			current_offset += stride_size * n_samples;
		}
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

}
}
