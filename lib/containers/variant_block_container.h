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
#include "core/variant_site_annotation.h"
#include "variant_block_mapper.h"
#include "core/variant_reader_objects.h"

namespace tachyon {
namespace containers {

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
    typedef core::VariantHeader   global_header_type;
	typedef containers::VariantBlock             block_entry_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef containers::GenotypeSummary          genotype_summary_type;
	typedef containers::VariantSiteAnnotation    site_annotation_type;
	typedef containers::IntervalContainer        interval_container_type;
	typedef HashContainer                        hash_container_type;
	typedef DataBlockSettings                    block_settings_type;
	typedef VariantReaderObjects                 objects_type;

public:
	VariantBlockContainer() :
		header_(nullptr)
	{

	}

	VariantBlockContainer(const global_header_type& header) :
		mapper_(header.header_magic.n_format_values, header.header_magic.n_info_values),
		header_(&header)
	{

	}

	~VariantBlockContainer(void){}

	self_type& operator<<(const global_header_type& header){
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

	bool loadObjects(block_settings_type& block_settings) const;

	/**<
	 * Primary construction function for generating the appropriate instances of
	 * iterators / containers
	 * @param objects Target objects
	 * @return        Returns reference to input target objects
	 */
	objects_type& loadObjects(objects_type& objects, block_settings_type& block_settings) const;

private:
	block_mapper_type    mapper_; // global -> local, local -> global, loaded or not, primitive type
	block_type           block_;
	hash_container_type  h_tables_format_;
	hash_container_type  h_tables_info_;
	site_annotation_type annotations_;
	const global_header_type* header_;
};

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
