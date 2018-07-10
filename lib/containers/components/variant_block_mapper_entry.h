#ifndef CONTAINERS_COMPONENTS_VARIANT_BLOCK_MAPPER_ENTRY_H_
#define CONTAINERS_COMPONENTS_VARIANT_BLOCK_MAPPER_ENTRY_H_

#include "data_container_header.h"

namespace tachyon{
namespace containers{

struct VariantBlockMapperEntry {
public:
	typedef VariantBlockMapperEntry         self_type;
	typedef containers::DataContainerHeader header_type;

public:
	VariantBlockMapperEntry() :
		load_order_index(-1),
		stream_id_local(-1),
		stream_id_global(-1),
		offset(nullptr)
	{}

	VariantBlockMapperEntry(const U32 load_order_index, const S32 target_stream_disk, const header_type* offset) :
		load_order_index(load_order_index),
		stream_id_local(target_stream_disk),
		stream_id_global(offset->data_header.global_key),
		offset(offset)
	{}

	~VariantBlockMapperEntry(){}

	VariantBlockMapperEntry(const VariantBlockMapperEntry& other) :
		load_order_index(other.load_order_index),
		stream_id_local(other.stream_id_local),
		stream_id_global(other.stream_id_global),
		offset(other.offset)
	{}

	VariantBlockMapperEntry(VariantBlockMapperEntry&& other) :
		load_order_index(other.load_order_index),
		stream_id_local(other.stream_id_local),
		stream_id_global(other.stream_id_global),
		offset(other.offset)
	{}

	VariantBlockMapperEntry& operator=(const VariantBlockMapperEntry& other){
		this->load_order_index = other.load_order_index;
		this->stream_id_local  = other.stream_id_local;
		this->stream_id_global = other.stream_id_global;
		this->offset = other.offset;
		return *this;
	}

	VariantBlockMapperEntry& operator=(VariantBlockMapperEntry&& other){
		if(this!=&other) // prevent self-move
		{
			this->load_order_index = other.load_order_index;
			this->stream_id_local  = other.stream_id_local;
			this->stream_id_global = other.stream_id_global;
			this->offset = other.offset;
		}
		return *this;
	}

	inline bool operator<(const self_type& other) const{ return(this->offset->data_header.offset < other.offset->data_header.offset); }
	inline bool operator>(const self_type& other) const{ return(!((*this) < other)); }

	void operator()(const U32& load_order_index, const U32& stream_id_local, const S32& stream_id_global, const header_type* offset){
		this->load_order_index = load_order_index;
		this->stream_id_local = stream_id_local;
		this->stream_id_global = stream_id_global;
		this->offset = offset;
	}

public:
	S32 load_order_index;      // Loaded order index
	S32 stream_id_local;       // Local target index
	S32 stream_id_global;      // Global target index
	const header_type* offset; // Header object of target data container
};

}
}



#endif /* CONTAINERS_COMPONENTS_VARIANT_BLOCK_MAPPER_ENTRY_H_ */
