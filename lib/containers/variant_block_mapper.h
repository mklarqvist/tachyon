#ifndef CONTAINERS_VARIANT_BLOCK_MAPPER_H_
#define CONTAINERS_VARIANT_BLOCK_MAPPER_H_

#include "components/variant_block_mapper_entry.h"
#include "components/variant_block_footer.h"
#include "core/header/variant_header.h"
#include "core/data_block_settings.h"

namespace tachyon{
namespace containers{

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
	inline bool empty(void) const{ return(this->values_.empty()); }
	inline size_type  size(void) const{ return(this->values_.size()); }

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


// Forward declaration
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
	typedef VariantHeader               header_type;
	typedef DataBlockSettings           block_settings_type;

public:
	VariantBlockMapper() ;
	VariantBlockMapper(const size_t n_format_fields, const size_t n_info_fields);
	~VariantBlockMapper();

	self_type& operator<<(const header_type& header){
		this->n_format_fields = header.format_fields_.size();
		this->n_info_fields   = header.info_fields_.size();
		return(*this);
	}

	/**<
	 * Update the current mapper object with the provided block footer
	 * @param block_footer        Target block footer providing information used for mappings
	 * @return                    Returns TRUE if successful or FALSE otherwise
	 */
	bool build(const block_footer_type& block_footer);

	/**<
	 * Update the current mapper object with the provided block footer, desired info keys
	 * and format keys
	 * @param info_keys    Provided info keys
	 * @param format_keys  Provided format keys
	 * @param block_footer Associated variant block footer information
	 * @return             Returns TRUE upon success or FALSE otherwise
	 */
	bool build(const std::vector<U32>& info_keys, const std::vector<U32>& format_keys, const block_footer_type& block_footer);

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
	inline size_t getNumberFormatLoaded(void) const{ return(this->format_container_loaded_.size()); }
	inline size_t getNumberInfoLoaded(void) const{ return(this->info_container_loaded_.size()); }

private:
	size_type n_format_fields;          // Total number of format fields in the file NOT the block
	size_type n_info_fields;            // Total number of info fields in the file NOT the block
	value_type info_container_global_;  // Global -> local mapping
	value_type info_container_local_;   // Local -> global mapping
	value_type info_container_loaded_;  // Loaded order -> local
	value_type format_container_global_;// Global -> local mapping
	value_type format_container_local_; // Local -> global mapping
	value_type format_container_loaded_;// Loaded order -> local
};

}
}



#endif /* CONTAINERS_VARIANT_BLOCK_MAPPER_H_ */
