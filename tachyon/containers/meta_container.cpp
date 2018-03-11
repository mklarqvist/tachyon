#include "meta_container.h"

namespace tachyon{
namespace containers{

MetaContainer::MetaContainer(const VariantBlock& block) :
	n_entries(block.meta_hot_container.buffer_data_uncompressed.size() / sizeof(hot_type)),
	__entries(nullptr)
{
	this->__ctor_setup(block);
}

MetaContainer::~MetaContainer(void){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__entries + i)->~MetaEntry();

	::operator delete[](static_cast<void*>(this->__entries));
}

void MetaContainer::__ctor_setup(const VariantBlock& block){
	// 1
	// cast hot_entries and determine n_entries
	assert(block.meta_hot_container.buffer_data_uncompressed.size() % sizeof(hot_type) == 0);
	const hot_type* const hot_entries = reinterpret_cast<const hot_type* const>(block.meta_hot_container.buffer_data_uncompressed.buffer);

	// allocate memory
	// iteratively call ctor for metaentry(cold,hot) or metaentry(hot)
	const U64& limit = block.meta_cold_container.buffer_data_uncompressed.size();
	if(limit){
		// Determine number of entries
		U32 current_offset = 0;
		U32 count_entries  = 0;

		while(true){
			const U32& l_body = *reinterpret_cast<const U32* const>(&block.meta_cold_container.buffer_data_uncompressed.buffer[current_offset]);
			assert(current_offset + l_body <= limit);
			current_offset += l_body;
			++count_entries;
			if(current_offset == limit) break;
		}
		assert(count_entries == this->n_entries);
		this->__entries = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

		current_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			const U32& l_body = *reinterpret_cast<const U32* const>(&block.meta_cold_container.buffer_data_uncompressed.buffer[current_offset]);
			new( &this->__entries[i] ) value_type( hot_entries[i], &block.meta_cold_container.buffer_data_uncompressed.buffer[current_offset] );
			current_offset += l_body;
		}
	} else {
		this->__entries = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
		for(U32 i = 0; i < this->n_entries; ++i){
			new( &this->__entries[i] ) value_type( hot_entries[i] );
		}
	}


	//U32 n_filter_ids = 0;
	//U32 n_format_ids = 0;
	//U32 n_info_ids   = 0;
	if(block.meta_filter_map_ids.buffer_data_uncompressed.size()){
		assert(block.meta_filter_map_ids.buffer_data_uncompressed.size() % sizeof(U32) == 0);
		const S32* const filter_ids = reinterpret_cast<const S32* const>(block.meta_filter_map_ids.buffer_data_uncompressed.data());
		if(block.meta_filter_map_ids.header.data_header.isUniform()){
			assert(block.meta_filter_map_ids.buffer_data_uncompressed.size() == sizeof(U32));
			//n_filter_ids = this->size();
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).filter_pattern_id = filter_ids[0];
		} else {
			assert(block.meta_filter_map_ids.buffer_data_uncompressed.size() / sizeof(U32) == this->size());
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).filter_pattern_id = filter_ids[i];

			//n_filter_ids = block.meta_filter_map_ids.buffer_data_uncompressed.size() / sizeof(U32);
		}
	}

	if(block.meta_info_map_ids.buffer_data_uncompressed.size()){
		assert(block.meta_info_map_ids.buffer_data_uncompressed.size() % sizeof(U32) == 0);
		const S32* const info_ids = reinterpret_cast<const S32* const>(block.meta_info_map_ids.buffer_data_uncompressed.data());
		if(block.meta_info_map_ids.header.data_header.isUniform()){
			assert(block.meta_info_map_ids.buffer_data_uncompressed.size() == sizeof(U32));
			//n_info_ids = this->size();
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).info_pattern_id = info_ids[0];

		} else {
			assert(block.meta_info_map_ids.buffer_data_uncompressed.size() / sizeof(U32) == this->size());
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).info_pattern_id = info_ids[i];

			//n_info_ids = block.meta_info_map_ids.buffer_data_uncompressed.size() / sizeof(U32);
		}
	}

	if(block.meta_format_map_ids.buffer_data_uncompressed.size()){
		assert(block.meta_format_map_ids.buffer_data_uncompressed.size() % sizeof(U32) == 0);
		const S32* const format_ids = reinterpret_cast<const S32* const>(block.meta_format_map_ids.buffer_data_uncompressed.data());
		if(block.meta_format_map_ids.header.data_header.isUniform()){
			assert(block.meta_format_map_ids.buffer_data_uncompressed.size() == sizeof(U32));
			//n_format_ids = this->size();
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).format_pattern_id = format_ids[0];

		} else {
			assert(block.meta_format_map_ids.buffer_data_uncompressed.size() / sizeof(U32) == this->size());
			for(U32 i = 0; i < this->size(); ++i)
				this->at(i).format_pattern_id = format_ids[i];

			//n_format_ids = block.meta_format_map_ids.buffer_data_uncompressed.size() / sizeof(U32);
		}
	}

	//assert(n_filter_ids == n_format_ids);
	//assert(n_format_ids == n_info_ids);
}

}
}
