#include "meta_cold_container.h"

namespace tachyon{
namespace containers{

MetaColdContainer::MetaColdContainer() : n_entries(0), __entries(nullptr){}

MetaColdContainer::MetaColdContainer(const VariantBlock& block) :
	n_entries(0),
	__entries(nullptr)
{
	this->__ctor_setup(block.meta_cold_container);
}

MetaColdContainer::MetaColdContainer(const DataContainer& container) :
	n_entries(0),
	__entries(nullptr)
{
	this->__ctor_setup(container);
}

MetaColdContainer::~MetaColdContainer(void){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__entries + i)->~MetaCold();

	::operator delete[](static_cast<void*>(this->__entries));
}

void MetaColdContainer::__ctor_setup(const DataContainer& container){
	// Determine number of entries
	U32 current_offset = 0;
	U32 count_entries  = 0;
	const U64& limit   = container.buffer_data_uncompressed.size();
	while(true){
		const U32& l_body = *reinterpret_cast<const U32* const>(&container.buffer_data_uncompressed.buffer[current_offset]);
		assert(current_offset + l_body <= limit);
		current_offset += l_body;
		++count_entries;
		if(current_offset == limit) break;
	}
	this->n_entries = count_entries;
	this->__entries = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	current_offset = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		const U32& l_body = *reinterpret_cast<const U32* const>(&container.buffer_data_uncompressed.buffer[current_offset]);
		new( &this->__entries[i] ) value_type( &container.buffer_data_uncompressed.buffer[current_offset] );
		current_offset += l_body;
	}
}

}
}
