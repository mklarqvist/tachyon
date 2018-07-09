#ifndef CONTAINERS_VARIANT_BLOCK_CONTAINER_H_
#define CONTAINERS_VARIANT_BLOCK_CONTAINER_H_

#include "variant_block.h"
#include "containers/meta_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/interval_container.h"
#include "containers/hash_container.h"

namespace tachyon {
namespace containers {

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



class VariantBlockMapperContainer {
private:
	typedef VariantBlockMapperContainer self_type;
	typedef VariantBlockMapperEntry     value_type;
	typedef std::size_t                 size_type;
	typedef value_type&                 reference;
	typedef const value_type&           const_reference;
	typedef value_type*                 pointer;
	typedef const value_type*           const_pointer;

public:
	VariantBlockMapperContainer(){}
	~VariantBlockMapperContainer(){}

    // Capacity
	inline const bool       empty(void) const{ return(this->values_.empty()); }
	inline const size_type& size(void) const{ return(this->size()); }

	// Element access
	inline pointer         data(void){ return(&this->values_[0]); }
	inline const_pointer   data(void) const{ return(&this->values_[0]); }
	inline reference       operator[](const U32& position){ return(this->values_[position]); }
	inline const_reference operator[](const U32& position) const{ return(this->values_[position]); }
	inline reference       at(const U32& position){ return(this->values_[position]); }
	inline const_reference at(const U32& position) const{ return(this->values_[position]); }

	inline void clear(void){ this->values_.clear(); }
	inline void resize(const size_t& new_size){ this->values_.resize(new_size); }

private:
	std::vector<value_type> values_;
};


/**<
 * Mapper class for VariantBlock: allows easy mapping from
 * Global -> local
 * Local  -> global
 * Loaded order -> local
 */
class VariantBlockMapper {
private:
	typedef VariantBlockMapper          self_type;
	typedef VariantBlockMapperContainer value_type;
	typedef std::size_t                 size_type;
	typedef value_type&                 reference;
	typedef const value_type&           const_reference;
	typedef value_type*                 pointer;
	typedef const value_type*           const_pointer;
	typedef VariantBlockMapperEntry     map_type;
	typedef VariantBlockFooter          block_footer_type;
	 typedef core::VariantHeader        header_type;

public:
	VariantBlockMapper() :
		n_format_fields(0),
		n_info_fields(0)
	{

	}

	VariantBlockMapper(const size_t n_format_fields, const size_t n_info_fields) :
		n_format_fields(n_format_fields),
		n_info_fields(n_info_fields)
	{

	}

	~VariantBlockMapper(){}

	self_type& operator<<(const header_type& header){
		this->n_format_fields = header.header_magic.n_format_values;
		this->n_info_fields   = header.header_magic.n_info_values;
		return(*this);
	}

	/**<
	 * Update the current mapper object with the provided block footer
	 * @param block_footer Target block footer used for mapping
	 * @return             Returns TRUE if successful or FALSE otherwise
	 */
	bool update(const block_footer_type& block_footer){
		this->format_container_global_.clear();
		this->format_container_global_.resize(this->n_format_fields);
		this->format_container_local_.clear();
		this->format_container_local_.resize(block_footer.n_format_streams);

		for(U32 i = 0; i < block_footer.n_format_streams; ++i){
			// Set global -> local mapping
			this->format_container_global_[block_footer.format_offsets[i].data_header.global_key](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
			// Set local -> global mapping
			this->format_container_local_[i](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
		}

		this->info_container_global_.clear();
		this->info_container_global_.resize(this->n_info_fields);
		this->info_container_local_.clear();
		this->info_container_local_.resize(block_footer.n_info_streams);
		for(U32 i = 0; i < block_footer.n_info_streams; ++i){
			// Set global -> local mapping
			this->info_container_global_[block_footer.info_offsets[i].data_header.global_key](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
			// Set local -> global mapping
			this->info_container_local_[i](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
		}

		return true;
	}

	inline bool isLoadedFormatGlobal(const U32& key) const;
	inline bool isLoadedInfoGlobal(const U32& key) const;
	inline bool isLoadedFormatLocal(const U32& key) const;
	inline bool isLoadedInfoLocal(const U32& key) const;
	map_type& getGlobalFormat(const U32& key);
	map_type& getLocalFormat(const U32& key);
	map_type& getGlobalInfo(const U32& key);
	map_type& getLocalInfo(const U32& key);

private:
	size_type n_format_fields;
	size_type n_info_fields;
	value_type info_container_global_;  // Global -> local mapping
	value_type info_container_local_;   // Local -> global mapping
	value_type info_container_loaded_;  // Loaded order -> local
	value_type format_container_global_;// Global -> local mapping
	value_type format_container_local_; // Local -> global mapping
	value_type format_container_loaded_;// Loaded order -> local
};


/**<
 * Decouples the low-level object `VariantBlock` that holds all internal data and
 * strides. This container should preferably be the primary object the end-user
 * interacts with.
 */
class VariantBlockContainer {
private:
	typedef VariantBlockContainer self_type;
    typedef DataContainerHeader   data_header_type;
    typedef VariantBlock          block_type;
    typedef VariantBlockHeader    block_header_type;
    typedef VariantBlockFooter    block_footer_type;
    typedef VariantBlockMapper    block_mapper_type;
    typedef core::VariantHeader   header_type;
	typedef containers::VariantBlock             block_entry_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef containers::GenotypeSummary          genotype_summary_type;
	typedef containers::IntervalContainer        interval_container_type;
	typedef HashContainer                        hash_container_type;

public:
	VariantBlockContainer()
	{

	}

	VariantBlockContainer(const header_type& header) :
		mapper_(header.header_magic.n_format_values, header.header_magic.n_info_values)
	{

	}

	~VariantBlockContainer(void){}

	self_type& operator<<(const header_type& header){
		this->mapper_ << header;
		return(*this);
	}

	void reset(void){
		this->block_.clear();
		// Todo: reset hashes
	}

	bool build(){
		if(this->mapper_.update(this->block_.footer) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to build mapper..." << std::endl;
			return false;
		}
		return(true);
	}

	inline block_type& getBlock(void){ return(this->block_); }
	inline const block_type& getBlock(void) const{ return(this->block_); }

	inline const bool anyEncrypted(void) const{ return(this->block_.header.controller.anyEncrypted); }

private:
	block_mapper_type    mapper_; // global -> local, local -> global, loaded or not, primitive type
	block_type           block_;
	hash_container_type  h_tables_format_;
	hash_container_type  h_tables_info_;
};

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
