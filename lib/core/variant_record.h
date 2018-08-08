#ifndef CORE_VARIANT_RECORD_H_
#define CORE_VARIANT_RECORD_H_

#include "containers/meta_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"

namespace tachyon{

// Forward declare.
namespace containers {
class VariantBlockContainer;
}

struct yon1_t {
	yon1_t(void) :
		is_dirty(false),
		is_loaded_meta(false),
		is_loaded_gt(false),
		n_format(0), n_info(0), n_filter(0),
		id_block(0),
		meta(nullptr),
		gt(nullptr),
		info(nullptr),
		fmt(nullptr),
		info_containers(nullptr),
		format_containers(nullptr),
		gt_i(nullptr),
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
		delete this->gt;
	}

	bool is_dirty; // if data has been modified in the raw buffer but not the containers
	bool is_loaded_meta;
	bool is_loaded_gt;
	uint16_t n_format, n_info, n_filter;
	uint32_t id_block; // incremental id in the block container
	core::MetaEntry* meta;
	yon_gt* gt;
	containers::PrimitiveContainerInterface** info;
	containers::PrimitiveGroupContainerInterface** fmt;
	containers::InfoContainerInterface** info_containers;
	containers::FormatContainerInterface** format_containers;
	containers::GenotypeContainerInterface* gt_i;
	std::vector<const YonInfo*> info_hdr;
	std::vector<const YonFormat*> format_hdr;
	std::vector<const io::VcfFilter*> filter_hdr;
	std::vector<int>* info_ids;
	std::vector<int>* format_ids;
	std::vector<int>* filter_ids;
	containers::VariantBlockContainer* parent_container;
};

}



#endif /* CORE_VARIANT_RECORD_H_ */
