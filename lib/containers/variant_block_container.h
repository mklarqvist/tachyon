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
#include "components/variant_block_mapper_entry.h"

namespace tachyon {
namespace containers {

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
	~VariantBlockMapperContainer() = default;

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

	// Operators
	inline void operator+=(const value_type& value){ this->values_.push_back(value); }

	// Utility
	inline void clear(void){ this->values_.clear(); }
	inline void resize(const size_t& new_size){ this->values_.resize(new_size); }

private:
	std::vector<value_type> values_;
};


class VariantBlockContainer;

/**<
 * Mapper class for VariantBlock: allows easy mapping from
 * Global -> local
 * Local  -> global
 * Loaded order -> local
 * Patterns -> keys
 */
class VariantBlockMapper {
	friend VariantBlockContainer;
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
	 * Update the current mapper object with the provided block footer, desired info keys
	 * and format keys
	 * @param info_keys    Provided info keys
	 * @param format_keys  Provided format keys
	 * @param block_footer Associated variant block footer information
	 * @return             Returns TRUE upon success or FALSE otherwise
	 */
	bool build(const std::vector<U32>& info_keys, const std::vector<U32>& format_keys, const block_footer_type& block_footer){
		this->format_container_global_.clear();
		this->format_container_global_.resize(this->n_format_fields);
		this->format_container_local_.clear();
		this->format_container_local_.resize(block_footer.n_format_streams);
		this->format_container_loaded_.clear();

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

	inline bool isLoadedFormatGlobal(const U32& key) const{ return(this->format_container_global_[key].load_order_index != -1); }
	inline bool isLoadedInfoGlobal(const U32& key) const{ return(this->info_container_global_[key].load_order_index != -1); }
	inline bool isLoadedFormatLocal(const U32& key) const{ return(this->format_container_local_[key].load_order_index != -1); }
	inline bool isLoadedInfoLocal(const U32& key) const{ return(this->info_container_local_[key].load_order_index != -1); }
	inline map_type& getGlobalFormat(const U32& key){ return(this->format_container_global_[key]); }
	inline map_type& getLocalFormat(const U32& key){ return(this->format_container_local_[key]); }
	inline map_type& getLoadedFormat(const U32& key){ return(this->format_container_loaded_[key]); }
	inline map_type& getGlobalInfo(const U32& key){ return(this->info_container_global_[key]); }
	inline map_type& getLocalInfo(const U32& key){ return(this->info_container_local_[key]); }
	inline map_type& getLoadedInfo(const U32& key){ return(this->info_container_loaded_[key]); }
	inline const map_type& getGlobalFormat(const U32& key) const{ return(this->format_container_global_[key]); }
	inline const map_type& getLocalFormat(const U32& key) const{ return(this->format_container_local_[key]); }
	inline const map_type& getLoadedFormat(const U32& key) const{ return(this->format_container_loaded_[key]); }
	inline const map_type& getGlobalInfo(const U32& key) const{ return(this->info_container_global_[key]); }
	inline const map_type& getLocalInfo(const U32& key) const{ return(this->info_container_local_[key]); }
	inline const map_type& getLoadedInfo(const U32& key) const{ return(this->info_container_loaded_[key]); }
	inline const size_t getNumberFormatLoaded(void) const{ return(this->format_container_loaded_.size()); }
	inline const size_t getNumberInfoLoaded(void) const{ return(this->info_container_loaded_.size()); }

private:
	size_type n_format_fields; // Total number of format fields in the file NOT the block
	size_type n_info_fields;   // Total number of info fields in the file NOT the block
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
	VariantBlockContainer() :
		header_(nullptr)
	{

	}

	VariantBlockContainer(const header_type& header) :
		mapper_(header.header_magic.n_format_values, header.header_magic.n_info_values),
		header_(&header)
	{

	}

	~VariantBlockContainer(void){}

	self_type& operator<<(const header_type& header){
		this->header_ = &header;
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
	bool readBlock(std::ifstream& stream, block_settings_type& settings);

	/**<
	 * Calculates which INFO pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name INFO field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_info_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Calculates which FORMAT pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name FORMAT field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_format_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Factory function for FORMAT container given an input `field` name
	 * @param field_name FORMAT field name to create a container for
	 * @return           Returns an instance of a `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_format_container(const std::string& field_name) const{
		core::HeaderMapEntry* match = nullptr;
		int format_field_global_id = -2;
		if(this->header_->getFormatField(field_name, match)){
			format_field_global_id = match->IDX;
		} else return nullptr;

		if(format_field_global_id >= 0){
			const S32& target_local_id = this->getMapper().getGlobalFormat(format_field_global_id).load_order_index;
			if(target_local_id < 0) return nullptr;
			return(new containers::FormatContainer<T>(this->getBlock().format_containers[target_local_id], this->header_->getSampleNumber()));
		}
		else return nullptr;
	}

	/**<
	 * Factory function for balanced FORMAT container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     FORMAT field name to create a container for
	 * @param meta_container Container for meta objects used to balance the FORMAT container
	 * @return               Returns an instance of a balanced `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_balanced_format_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		if(meta_container.size() == 0)
			return new containers::FormatContainer<T>();

		core::HeaderMapEntry* match = nullptr;
		int format_field_global_id = -2;
		if(this->header_->getFormatField(field_name, match)){
			format_field_global_id = match->IDX;
		} else return nullptr;

		if(format_field_global_id >= 0){
			const std::vector<bool> pattern_matches = this->get_format_field_pattern_matches(field_name);
			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			const S32& target_local_id = this->getMapper().getGlobalFormat(format_field_global_id).load_order_index;
			if(target_local_id < 0) return nullptr;

			return(new containers::FormatContainer<T>(this->getBlock().format_containers[target_local_id], meta_container, pattern_matches, this->header_->getSampleNumber()));
		}
		else return nullptr;
	}

	/**<
	 * Factory function for INFO container given an input `field` name
	 * @param field_name INFO field name to create a container for
	 * @return           Returns an instance of a `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_info_container(const std::string& field_name) const{
		core::HeaderMapEntry* match = nullptr;
		int info_field_global_id = -2;
		if(this->header_->getInfoField(field_name, match)){
			info_field_global_id = match->IDX;
		} else return nullptr;

		if(info_field_global_id >= 0){
			const S32& target_local_id = this->getMapper().getGlobalInfo(info_field_global_id).load_order_index;
			if(target_local_id < 0) return nullptr;
			return(new containers::InfoContainer<T>(this->getBlock().info_containers[target_local_id]));
		}
		else return nullptr;
	}

	/**<
	 * Factory function for balanced INFO container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     INFO field name to create a container for
	 * @param meta_container Container for meta objects used to balance the INFO container
	 * @return               Returns an instance of a balanced `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_balanced_info_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		if(meta_container.size() == 0)
			return new containers::InfoContainer<T>();

		core::HeaderMapEntry* match = nullptr;
		int info_field_global_id = -2;
		if(this->header_->getInfoField(field_name, match)){
			info_field_global_id = match->IDX;
		} else return nullptr;

		if(info_field_global_id >= 0){
			const std::vector<bool> pattern_matches = this->get_info_field_pattern_matches(field_name);

			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			const S32& target_local_id = this->getMapper().getGlobalInfo(info_field_global_id).load_order_index;
			if(target_local_id < 0) return nullptr;

			return(new containers::InfoContainer<T>(this->getBlock().info_containers[target_local_id], meta_container, pattern_matches));
		}
		else return nullptr;
	}

	// Accessors
	inline block_type& getBlock(void){ return(this->block_); }
	inline const block_type& getBlock(void) const{ return(this->block_); }
	inline block_mapper_type& getMapper(void){ return(this->mapper_); }
	inline const block_mapper_type& getMapper(void) const{ return(this->mapper_); }

	// Checkers
	inline const bool anyEncrypted(void) const{ return(this->block_.header.controller.anyEncrypted); }
	inline const bool hasGenotypes(void) const{ return(this->block_.header.controller.hasGT); }
	inline const bool hasPermutedGenotypes(void) const{ return(this->block_.header.controller.hasGTPermuted); }

private:
	block_mapper_type    mapper_; // global -> local, local -> global, loaded or not, primitive type
	block_type           block_;
	hash_container_type  h_tables_format_;
	hash_container_type  h_tables_info_;
	const header_type*   header_;
};

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
