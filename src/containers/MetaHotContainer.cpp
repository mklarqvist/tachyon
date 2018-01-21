#include "MetaHotContainer.h"

namespace Tachyon{
namespace Core{

MetaHotContainer::MetaHotContainer() : n_entries(0), __entries(nullptr){}

MetaHotContainer::MetaHotContainer(const Block& block) :
	n_entries(block.meta_hot_container.buffer_data_uncompressed.size() / sizeof(value_type)),
	__entries(new value_type[this->n_entries])
{
	assert(block.meta_hot_container.buffer_data_uncompressed.size() % sizeof(value_type) == 0);

	const MetaHot* const d = reinterpret_cast<const MetaHot* const>(block.meta_hot_container.buffer_data_uncompressed.data);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->__entries[i] = d[i];
}

MetaHotContainer::MetaHotContainer(const Container& container) :
	n_entries(container.buffer_data_uncompressed.size() / sizeof(value_type)),
	__entries(new value_type[this->n_entries])
{
	assert(container.buffer_data_uncompressed.size() % sizeof(value_type) == 0);

	const MetaHot* const d = reinterpret_cast<const MetaHot* const>(container.buffer_data_uncompressed.data);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->__entries[i] = d[i];
}

MetaHotContainer::~MetaHotContainer(void){
	delete [] this->__entries;
}

}
}
