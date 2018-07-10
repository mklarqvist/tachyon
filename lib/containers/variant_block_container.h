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
	inline const size_type  size(void) const{ return(this->values_.size()); }

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
	typedef core::VariantHeader         header_type;
	typedef DataBlockSettings           block_settings_type;

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
	 * @param block_footer        Target block footer providing information used for mappings
	 * @return                    Returns TRUE if successful or FALSE otherwise
	 */
	bool build(const block_footer_type& block_footer){
		this->format_container_global_.clear();
		this->format_container_global_.resize(this->n_format_fields);
		this->format_container_local_.clear();
		this->format_container_local_.resize(block_footer.n_format_streams);

		for(U32 i = 0; i < block_footer.n_format_streams; ++i){
			// Set global -> local mapping -> loaded mapping
			this->format_container_global_[block_footer.format_offsets[i].data_header.global_key](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
			// Set local -> global mapping -> loaded mapping
			this->format_container_local_[i](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
		}

		this->info_container_global_.clear();
		this->info_container_global_.resize(this->n_info_fields);
		this->info_container_local_.clear();
		this->info_container_local_.resize(block_footer.n_info_streams);

		for(U32 i = 0; i < block_footer.n_info_streams; ++i){
			// Set global -> local mapping -> loaded mapping
			this->info_container_global_[block_footer.info_offsets[i].data_header.global_key](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
			// Set local -> global mapping -> loaded mapping
			this->info_container_local_[i](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
		}


		return true;
	}

	/**<
	 *
	 * @param info_keys
	 * @param format_keys
	 * @param block_footer
	 * @return
	 */
	bool build(const std::vector<U32>& info_keys, const std::vector<U32>& format_keys, const block_footer_type& block_footer){
		this->format_container_global_.clear();
		this->format_container_global_.resize(this->n_format_fields);
		this->format_container_local_.clear();
		this->format_container_local_.resize(block_footer.n_format_streams);
		this->format_container_loaded_.clear();

		// Map local first
		for(U32 i = 0; i < block_footer.n_format_streams; ++i){
			// Set global -> local mapping -> loaded mapping
			this->format_container_global_[block_footer.format_offsets[i].data_header.global_key](-1, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
			// Set local -> global mapping -> loaded mapping
			this->format_container_local_[i](-1, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
		}

		for(U32 i = 0; i < format_keys.size(); ++i){
			this->format_container_global_[format_keys[i]].load_order_index = i;
			this->format_container_local_[this->format_container_global_[format_keys[i]].stream_id_local].load_order_index = i;
		}

		this->info_container_global_.clear();
		this->info_container_global_.resize(this->n_info_fields);
		this->info_container_local_.clear();
		this->info_container_local_.resize(block_footer.n_info_streams);
		this->info_container_loaded_.clear();

		for(U32 i = 0; i < block_footer.n_info_streams; ++i){
			// Set global -> local mapping -> loaded mapping
			this->info_container_global_[block_footer.info_offsets[i].data_header.global_key](-1, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
			// Set local -> global mapping -> loaded mapping
			this->info_container_local_[i](-1, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
		}

		for(U32 i = 0; i < info_keys.size(); ++i){
			this->info_container_global_[info_keys[i]].load_order_index = i;
			this->info_container_local_[this->info_container_global_[info_keys[i]].stream_id_local].load_order_index = i;
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

public:
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
	typedef DataBlockSettings                    block_settings_type;

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

	bool buildMapper(const std::vector<U32>& info_keys, const std::vector<U32>& format_keys, const block_settings_type& block_settings){
		if(this->mapper_.build(info_keys, format_keys, this->block_.footer) == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to build mapper..." << std::endl;
			return false;
		}
		return(true);
	}

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity until later.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool readBlock(std::ifstream& stream, block_settings_type& settings){
		// Get info and format keys
		std::vector<U32> info_keys, format_keys;
		if(settings.info_ID_list.size()) info_keys = this->block_.intersectInfoKeys(settings.info_ID_list);
		else info_keys = this->block_.getInfoKeys();
		if(settings.format_ID_list.size()) format_keys = this->block_.intersectFormatKeys(settings.format_ID_list);
		else format_keys = this->block_.getFormatKeys();

		if(this->buildMapper(info_keys, format_keys, settings) == false)
			return false;

		// Todo: ascertain random access order is guaranteed

		if(settings.ppa.load){
			if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
				this->block_.ppa_manager.header = this->block_.footer.offset_ppa;
				stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offset_ppa.data_header.offset);
				stream >> this->block_.ppa_manager;
			}
		}

		if(settings.contig.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_contig, this->block_.meta_contig_container);
		}

		if(settings.positions.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_position, this->block_.meta_positions_container);
		}

		if(settings.controller.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_controllers, this->block_.meta_controller_container);
		}

		if(settings.quality.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_quality, this->block_.meta_quality_container);
		}

		if(settings.names.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_names, this->block_.meta_names_container);
		}

		if(settings.alleles.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_refalt, this->block_.meta_refalt_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_meta_alleles, this->block_.meta_alleles_container);
		}

		if(settings.genotypes_rle.load || settings.genotypes_all.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_8b, this->block_.gt_rle8_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_16b, this->block_.gt_rle16_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_32b, this->block_.gt_rle32_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_64b, this->block_.gt_rle64_container);
		}

		if(settings.genotypes_simple.load || settings.genotypes_all.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_simple8, this->block_.gt_simple8_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple16, this->block_.gt_simple16_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple32, this->block_.gt_simple32_container);
			this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple64, this->block_.gt_simple64_container);
		}

		if(settings.genotypes_support.load || settings.genotypes_all.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_helper, this->block_.gt_support_data_container);
		}

		if(settings.set_membership.load || settings.genotypes_all.load){
			this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_info_id, this->block_.meta_info_map_ids);
			this->block_.__loadContainer(stream, this->block_.footer.offset_meta_filter_id, this->block_.meta_filter_map_ids);
			this->block_.__loadContainer(stream, this->block_.footer.offset_meta_format_id, this->block_.meta_format_map_ids);
		}

		// Load all info
		if(settings.info_all.load && this->block_.footer.n_info_streams){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

			this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
			for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
				this->block_.__loadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
				this->mapper_.info_container_loaded_.at(i)(i, i, this->block_.footer.info_offsets[i].data_header.global_key, &this->block_.footer.info_offsets[i]);
			}
		}
		// If we have supplied a list of identifiers
		else if(settings.info_ID_list.size()){
			this->mapper_.info_container_loaded_.resize(info_keys.size());
			// Ascertain that random access is linearly forward
			for(U32 i = 0; i < info_keys.size(); ++i){
				stream.seekg(this->block_.start_compressed_data_ + this->mapper_.info_container_global_[info_keys[i]].offset->data_header.offset);
				if(!stream.good()){
					std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for INFO field in block!" << std::endl;
					return false;
				}

				this->block_.__loadContainer(stream, this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local], this->block_.info_containers[i]);
				this->mapper_.info_container_loaded_.at(i)(i, this->mapper_.info_container_global_[info_keys[i]].stream_id_local, info_keys[i], &this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local]);
			}
		} // end case load_info_ID

		// Load all FORMAT data
		if(settings.format_all.load && this->block_.footer.n_format_streams){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
			this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
			for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
				this->block_.__loadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
				this->mapper_.format_container_loaded_.at(i)(i, i, this->block_.footer.format_offsets[i].data_header.global_key, &this->block_.footer.format_offsets[i]);
			}
			//std::cerr << this->end_compressed_data_ << '/' << (U64)stream.tellg() << std::endl;
			assert(this->end_compressed_data_ == (U64)stream.tellg());
		} // If we have supplied a list of identifiers
		else if(settings.format_ID_list.size()){
			this->mapper_.format_container_loaded_.resize(format_keys.size());
			// Ascertain that random access is linearly forward
			for(U32 i = 0; i < format_keys.size(); ++i){
				stream.seekg(this->block_.start_compressed_data_ + this->mapper_.format_container_global_[format_keys[i]].offset->data_header.offset);
				if(!stream.good()){
					std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for FORMAT field in block!" << std::endl;
					return false;
				}

				this->block_.__loadContainer(stream, this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local], this->block_.format_containers[i]);
				this->mapper_.format_container_loaded_.at(i)(i, this->mapper_.format_container_global_[format_keys[i]].stream_id_local, format_keys[i], &this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local]);
			}
		} // end case load_info_ID

		stream.seekg(this->block_.end_block_); // seek to end-of-block
		return(true);
	}

	// Accessors
	inline block_type& getBlock(void){ return(this->block_); }
	inline const block_type& getBlock(void) const{ return(this->block_); }

	// Checkers
	inline const bool anyEncrypted(void) const{ return(this->block_.header.controller.anyEncrypted); }
	inline const bool hasGenotypes(void) const{ return(this->block_.header.controller.hasGT); }
	inline const bool hasPermutedGenotypes(void) const{ return(this->block_.header.controller.hasGTPermuted); }

public:
	block_mapper_type    mapper_; // global -> local, local -> global, loaded or not, primitive type
	block_type           block_;
	hash_container_type  h_tables_format_;
	hash_container_type  h_tables_info_;
};

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
