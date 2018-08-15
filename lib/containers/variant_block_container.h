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
#include "core/variant_reader_objects.h"
#include "core/variant_record.h"

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
    typedef VariantHeader         global_header_type;
	typedef containers::VariantBlock             block_entry_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef containers::IntervalContainer        interval_container_type;
	typedef DataBlockSettings                    block_settings_type;
	typedef VariantReaderObjects                 objects_type;

	typedef std::unordered_map<int, int>    map_type;

public:
	VariantBlockContainer() :
		loaded_genotypes(false),
		header_(nullptr),
		gt_exp(nullptr)
	{

	}

	VariantBlockContainer(const global_header_type& header) :
		loaded_genotypes(false),
		header_(&header),
		gt_exp(nullptr)
	{
		delete [] this->gt_exp;
	}

	~VariantBlockContainer(void){
	}

	self_type& operator<<(const global_header_type& header){
		this->header_ = &header;
		return(*this);
	}

	void reset(void){
		this->block_.clear();
	}

	bool ParseSettings(block_settings_type& settings);
	bool ParseLoadedPatterns(block_settings_type& settings);

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity until later.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool ReadBlock(std::ifstream& stream, block_settings_type& settings);

	/**<
	 * Factory function for FORMAT container given an input `field` name
	 * @param field_name FORMAT field name to create a container for
	 * @return           Returns an instance of a `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_format_container(const std::string& field_name) const;

	/**<
	 * Factory function for balanced FORMAT container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     FORMAT field name to create a container for
	 * @param meta_container Container for meta objects used to balance the FORMAT container
	 * @return               Returns an instance of a balanced `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_balanced_format_container(const std::string& field_name, const containers::MetaContainer& meta_container) const;

	/**<
	 * Factory function for INFO container given an input `field` name
	 * @param field_name INFO field name to create a container for
	 * @return           Returns an instance of a `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_info_container(const std::string& field_name) const;

	/**<
	 * Factory function for balanced INFO container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     INFO field name to create a container for
	 * @param meta_container Container for meta objects used to balance the INFO container
	 * @return               Returns an instance of a balanced `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_balanced_info_container(const std::string& field_name, const containers::MetaContainer& meta_container) const;

	// Accessors
	inline block_type& GetBlock(void){ return(this->block_); }
	inline const block_type& GetBlock(void) const{ return(this->block_); }

	// Checkers
	inline bool AnyEncrypted(void) const{ return(this->block_.header.controller.anyEncrypted); }
	inline bool HasGenotypes(void) const{ return(this->block_.header.controller.hasGT); }
	inline bool HasPermutedGenotypes(void) const{ return(this->block_.header.controller.hasGTPermuted); }

	inline void AllocateGenotypeMemory(void){ this->gt_exp = new yon_gt_rcd*[this->header_->GetNumberSamples()]; }
	inline yon_gt_rcd** GetAllocatedGenotypeMemory(void){ return(this->gt_exp); }
	inline yon_gt_rcd** GetAllocatedGenotypeMemory(void) const{ return(this->gt_exp); }

	/**<
	 * Primary construction function for generating the appropriate instances of
	 * iterators / containers
	 * @param objects Target objects
	 * @return        Returns reference to input target objects
	 */
	objects_type* LoadObjects(objects_type* objects, block_settings_type& block_settings);

	objects_type* LoadObjects(block_settings_type& block_settings){
		objects_type* obj = new objects_type();
		return(this->LoadObjects(obj, block_settings));
	}

	yon1_t* LazyEvaluate(objects_type& objects);

private:
	bool loaded_genotypes;
	block_type                block_;
	const global_header_type* header_;
	std::vector<int>          info_id_local_loaded;
	std::vector<int>          format_id_local_loaded;
	std::vector<int>          info_id_global_loaded;
	std::vector<int>          format_id_global_loaded;
	std::vector< std::vector<int> > info_patterns_local;
	std::vector< std::vector<int> > format_patterns_local;
	map_type info_map_global;
	map_type format_map_global;
	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	yon_gt_rcd** gt_exp;
};


// IMPLEMENTATION -------------------------------------------------------------



template <class T>
containers::FormatContainer<T>* VariantBlockContainer::get_format_container(const std::string& field_name) const{
	int format_field_global_id = -2;
	YonFormat* fmt = this->header_->GetFormat(field_name);
	if(fmt != nullptr){
		format_field_global_id = fmt->idx;
	} else return new containers::FormatContainer<T>();

	DataContainer* container = this->block_.GetFormatContainer(format_field_global_id);
	if(container == nullptr) new containers::FormatContainer<T>();

	return(new containers::FormatContainer<T>(*container, this->header_->GetNumberSamples()));
}

template <class T>
containers::FormatContainer<T>* VariantBlockContainer::get_balanced_format_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
	if(meta_container.size() == 0)
		return new containers::FormatContainer<T>();

	int format_field_global_id = -2;
	YonFormat* fmt = this->header_->GetFormat(field_name);
	if(fmt != nullptr){
		format_field_global_id = fmt->idx;
	} else return new containers::FormatContainer<T>();

	DataContainer* container = this->block_.GetFormatContainer(format_field_global_id);
	if(container == nullptr) new containers::FormatContainer<T>();

	const std::vector<bool> pattern_matches = this->block_.FormatPatternSetMembership(format_field_global_id);

	U32 matches = 0;
	for(U32 i = 0; i < pattern_matches.size(); ++i)
		matches += pattern_matches[i];

	if(matches == 0) return new containers::FormatContainer<T>();

	return(new containers::FormatContainer<T>(*container, meta_container, pattern_matches, this->header_->GetNumberSamples()));

}

template <class T>
containers::InfoContainer<T>* VariantBlockContainer::get_info_container(const std::string& field_name) const{
	int info_field_global_id = -2;
	YonInfo* info = this->header_->GetInfo(field_name);
	if(info != nullptr){
		info_field_global_id = info->idx;
	} else return new containers::InfoContainer<T>();

	DataContainer* container = this->block_.GetInfoContainer(info_field_global_id);
	if(container == nullptr) new containers::InfoContainer<T>();

	return(new containers::InfoContainer<T>(*container));
}

template <class T>
containers::InfoContainer<T>* VariantBlockContainer::get_balanced_info_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
	int info_field_global_id = -2;
	YonInfo* info = this->header_->GetInfo(field_name);
	if(info != nullptr){
		info_field_global_id = info->idx;
	} else return new containers::InfoContainer<T>();

	DataContainer* container = this->block_.GetInfoContainer(info_field_global_id);
	if(container == nullptr) new containers::InfoContainer<T>();

	const std::vector<bool> pattern_matches = this->block_.InfoPatternSetMembership(info_field_global_id);

	U32 matches = 0;
	for(U32 i = 0; i < pattern_matches.size(); ++i)
		matches += pattern_matches[i];

	if(matches == 0) return new containers::InfoContainer<T>();

	return(new containers::InfoContainer<T>(*container, meta_container, pattern_matches));
}

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
