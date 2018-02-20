#include "checksum_container.h"

namespace tachyon{
namespace containers{

ChecksumContainer::ChecksumContainer(void) :
		n_entries(0),
		n_capacity(1000),
		__entries(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::ChecksumContainer(const size_type capacity) :
		n_entries(0),
		n_capacity(capacity),
		__entries(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::~ChecksumContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__entries + i)->~DigitalDigestPair)();

	::operator delete[](static_cast<void*>(this->__entries));
}

bool ChecksumContainer::allocate(const size_type n_entries){
	if(n_entries > this->capacity())
		return false;

	this->n_entries = n_entries;
	for(size_type i = 0; i < this->size(); ++i){
		new( &this->__entries[this->n_entries] ) value_type( );
		this->__entries[i].compressed.initialize();
		this->__entries[i].uncompressed.initialize();
	}

	return true;
}

void ChecksumContainer::operator+=(const_reference value){
	if(this->size() + 1 == this->capacity()){
		// is full
		std::cerr << "is full: need resize" << std::endl;
		return;
	}

	// Invoke copy-ctor
	new( &this->__entries[this->n_entries] ) value_type( value );
	this->__entries[this->n_entries].compressed.initialize();
	this->__entries[this->n_entries].uncompressed.initialize();
	++this->n_entries;
}

void ChecksumContainer::finalize(void){
	if(this->size() == 0) return;
	for(size_type i = 0; i < this->size(); ++i)
		this->at(i).finalize();
}

bool ChecksumContainer::update(const block_type& block, const U32* const mapTable){
	for(U32 i = 0; i < block.index_entry.n_info_streams; ++i){
		assert(mapTable[block.index_entry.info_offsets[i].key] < this->size());
		if(!(*this)[mapTable[block.index_entry.info_offsets[i].key]].uncompressed.update(block.info_containers[i].buffer_data_uncompressed, block.info_containers[i].buffer_strides_uncompressed, block.info_containers[i].header.hasMixedStride())){
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}

		assert(block.index_entry.info_offsets[i].key < this->size());
		if(!(*this)[mapTable[block.index_entry.info_offsets[i].key]].compressed.update(block.info_containers[i].buffer_data, block.info_containers[i].buffer_strides, block.info_containers[i].header.hasMixedStride())){
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}
	}

	for(U32 i = 0; i < block.index_entry.n_format_streams; ++i){
		if(!(*this)[mapTable[block.index_entry.format_offsets[i].key]].uncompressed.update(block.format_containers[i].buffer_data_uncompressed, block.format_containers[i].buffer_strides_uncompressed, block.format_containers[i].header.hasMixedStride())){
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}

		if(!(*this)[mapTable[block.index_entry.format_offsets[i].key]].compressed.update(block.format_containers[i].buffer_data, block.format_containers[i].buffer_strides, block.format_containers[i].header.hasMixedStride())){
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}
	}
	return true;
}

}

}
