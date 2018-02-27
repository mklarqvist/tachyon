#include "meta_hot_container.h"

namespace tachyon{
namespace containers{

MetaHotContainer::MetaHotContainer() : n_entries(0), __entries(nullptr){}

MetaHotContainer::MetaHotContainer(const VariantBlock& block) :
	n_entries(block.meta_hot_container.buffer_data_uncompressed.size() / sizeof(value_type)),
	__entries(new value_type[this->n_entries])
{
	assert(block.meta_hot_container.buffer_data_uncompressed.size() % sizeof(value_type) == 0);

	const value_type* const d = reinterpret_cast<const value_type* const>(block.meta_hot_container.buffer_data_uncompressed.buffer);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->__entries[i] = d[i];
}

MetaHotContainer::MetaHotContainer(const DataContainer& container) :
	n_entries(container.buffer_data_uncompressed.size() / sizeof(value_type)),
	__entries(new value_type[this->n_entries])
{
	assert(container.buffer_data_uncompressed.size() % sizeof(value_type) == 0);

	const value_type* const d = reinterpret_cast<const value_type* const>(container.buffer_data_uncompressed.buffer);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->__entries[i] = d[i];
}

MetaHotContainer::~MetaHotContainer(void){
	delete [] this->__entries;
}

}
}
