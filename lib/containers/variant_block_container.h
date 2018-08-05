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
#include "core/variant_site_annotation.h"
#include "core/variant_reader_objects.h"

namespace tachyon {
namespace containers {

// Forward declare.
class VariantBlockContainer;

// Todo:
struct yon1_t {
	yon1_t(void) :
		is_dirty(false),
		is_loaded_meta(false),
		is_loaded_gt(false),
		id_block(0),
		meta(nullptr),
		info(nullptr),
		fmt(nullptr),
		info_containers(nullptr),
		format_containers(nullptr),
		gt(nullptr),
		info_ids(nullptr),
		format_ids(nullptr),
		filter_ids(nullptr),
		parent_container(nullptr)
	{

	}
	~yon1_t(void){
		delete [] this->info;
		delete [] this->fmt;
		delete [] this->info_containers;
		delete [] this->format_containers;
	}

	bool is_dirty; // if data has been modified in the raw buffer but not the containers
	bool is_loaded_meta;
	bool is_loaded_gt;
	uint32_t id_block; // incremental id in the block container
	core::MetaEntry* meta;
	PrimitiveContainerInterface** info;
	PrimitiveGroupContainerInterface** fmt;
	InfoContainerInterface** info_containers;
	FormatContainerInterface** format_containers;
	GenotypeContainerInterface* gt;
	std::vector<const YonInfo*> info_hdr;
	std::vector<const YonFormat*> format_hdr;
	std::vector<const io::VcfFilter*> filter_hdr;
	std::vector<int>* info_ids;
	std::vector<int>* format_ids;
	std::vector<int>* filter_ids;
	VariantBlockContainer* parent_container;
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
    typedef VariantHeader         global_header_type;
	typedef containers::VariantBlock             block_entry_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef containers::GenotypeSummary          genotype_summary_type;
	typedef containers::VariantSiteAnnotation    site_annotation_type;
	typedef containers::IntervalContainer        interval_container_type;
	typedef DataBlockSettings                    block_settings_type;
	typedef VariantReaderObjects                 objects_type;

public:
	VariantBlockContainer() :
		header_(nullptr),
		objects_(nullptr)
	{

	}

	VariantBlockContainer(const global_header_type& header) :
		header_(&header),
		objects_(nullptr)
	{

	}

	~VariantBlockContainer(void){
		delete this->objects_;
	}

	self_type& operator<<(const global_header_type& header){
		this->header_ = &header;
		return(*this);
	}

	void reset(void){
		this->block_.clear();
		delete this->objects_;
		this->objects_ = nullptr;
		// Todo: reset hashes
	}

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
	inline objects_type* GetObjects(void){ return(this->objects_); }

	// Checkers
	inline bool AnyEncrypted(void) const{ return(this->block_.header.controller.anyEncrypted); }
	inline bool HasGenotypes(void) const{ return(this->block_.header.controller.hasGT); }
	inline bool HasPermutedGenotypes(void) const{ return(this->block_.header.controller.hasGTPermuted); }

	/**<
	 * Primary construction function for generating the appropriate instances of
	 * iterators / containers
	 * @param objects Target objects
	 * @return        Returns reference to input target objects
	 */
	objects_type& LoadObjects(objects_type& objects, block_settings_type& block_settings);

	objects_type* LoadObjects(block_settings_type& block_settings){
		delete this->objects_;
		this->objects_ = new objects_type;
		this->LoadObjects(*this->objects_, block_settings);
		return this->objects_;
	}

	yon1_t* LazyEvaluate(objects_type& objects){
		// Lazy evaluation of interleaved data.
		yon1_t* records = new yon1_t[objects.meta_container->size()];

		// Iterate over the sites described in the meta container.
		for(U32 i = 0; i < objects.meta_container->size(); ++i){
			records[i].is_dirty = false;
			records[i].is_loaded_meta = true;
			records[i].is_loaded_gt = true;
			records[i].id_block = i;
			records[i].meta = &objects.meta_container->at(i);

			if(records[i].is_loaded_gt)
				records[i].gt = &objects.genotype_container->at(i);

			records[i].info_ids = &this->block_.footer.info_patterns[objects.meta_container->at(i).info_pattern_id].pattern;
			records[i].info = new PrimitiveContainerInterface*[objects.n_loaded_info];
			records[i].info_containers = new InfoContainerInterface*[objects.n_loaded_info];

			// Populate Info data.
			for(U32 j = 0; j < records[i].info_ids->size(); ++j){
				InfoContainerInterface* info_cnt = objects.info_container_map[this->header_->info_fields_[records[i].info_ids->at(j)].id];
				records[i].info_containers[j] = info_cnt;
				records[i].info_hdr.push_back(&this->header_->info_fields_[records[i].info_ids->at(j)]);

				//std::cerr << "INFO " << j << "/" << records[i].info_ids->size() << ": target: " << this->header_->info_fields_[records[i].info_ids->at(j)].id << " container entries = " << info_cnt->size() << std::endl;

				switch(this->header_->info_fields_[records[i].info_ids->at(j)].yon_type){
				case(YON_VCF_HEADER_INTEGER):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<S32>*>(info_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_FLOAT):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<float>*>(info_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_STRING):
				case(YON_VCF_HEADER_CHARACTER):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<std::string>*>(info_cnt)->at(i);
					break;
				default:
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<U32>*>(info_cnt)->at(0);
					break;
				}
				//std::cerr << "target entries :" << records[i].info[j]->size() << std::endl;
			}

			//std::cerr << "start format" << std::endl;
			records[i].format_ids = &this->block_.footer.format_patterns[objects.meta_container->at(i).format_pattern_id].pattern;
			records[i].fmt = new PrimitiveGroupContainerInterface*[objects.n_loaded_info];
			records[i].format_containers = new FormatContainerInterface*[objects.n_loaded_format];

			//std::cerr << "after base format" << std::endl;

			// Populate Format data.
			for(U32 j = 0; j < records[i].format_ids->size(); ++j){
				FormatContainerInterface* fmt_cnt = objects.format_container_map[this->header_->format_fields_[records[i].format_ids->at(j)].id];
				records[i].format_containers[j] = fmt_cnt;
				records[i].format_hdr.push_back(&this->header_->format_fields_[records[i].format_ids->at(j)]);


				//std::cerr << "FORMAT " << j << "/" << records[i].format_ids->size() << ": target: " << this->header_->format_fields_[records[i].format_ids->at(j)].id << " container entries = " << fmt_cnt->size() << std::endl;

				// Todo: if format container is completely empty then skip
				if(fmt_cnt->size() == 0){
					records[i].fmt[j] = nullptr;
					continue;
				}

				switch(this->header_->format_fields_[records[i].format_ids->at(j)].yon_type){
				case(YON_VCF_HEADER_INTEGER):
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<S32>*>(fmt_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_FLOAT):
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<float>*>(fmt_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_STRING):
				case(YON_VCF_HEADER_CHARACTER):
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<std::string>*>(fmt_cnt)->at(i);
					break;
				default:
					std::cerr << "fmt at default" << std::endl;
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<U32>*>(fmt_cnt)->at(0);
					break;
				}
				//std::cerr << "target entries :" << records[i].fmt[j]->size() << std::endl;
			}

			// Populate Filter.
			records[i].filter_ids = &this->block_.footer.filter_patterns[objects.meta_container->at(i).filter_pattern_id].pattern;
			for(U32 j = 0; j < records[i].filter_ids->size(); ++j){
				records[i].filter_hdr.push_back(&this->header_->filter_fields_[records[i].filter_ids->at(j)]);
			}

			records[i].parent_container = this;
		}

		//std::cerr << "done" << std::endl;
		return(records);
	}

private:
	block_type                block_;
	site_annotation_type      annotations_;
	const global_header_type* header_;
	objects_type*             objects_;
	std::vector<int>          info_id_loaded;
	std::vector<int>          format_id_loaded;
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
