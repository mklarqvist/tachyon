#ifndef CORE_VARIANT_RECORD_H_
#define CORE_VARIANT_RECORD_H_

#include "containers/meta_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "occ.h"

namespace tachyon{

// Forward declare to allow variant to reference
// parent (host) container.
namespace containers {
class VariantBlockContainer;
}

/**<
 * Primary lazy-evaluated record of a variant. Construction is done
 * outside of this definition. The struct is constructed from pointers
 * to the appropriate location in other lazy-evaluated objects such as
 * Info, Format, Genotype, and Meta containers. Evaluation into this object
 * is relatively inexpensive in isolation but quiet expensive when considered
 * across the an entire dataset as this structure needs be evaluated many
 * millions of times over many different containers.
 */
struct yon1_t {
public:
	yon1_t(void);
	~yon1_t(void);
	bool EvaluateSummary(bool lazy_evaluate = true);
	bool EvaluateOcc();

public:
	bool is_dirty; // if data has been modified in the raw buffer but not the containers
	bool is_loaded_meta;
	bool is_loaded_gt;
	uint16_t n_format, n_info, n_filter;
	uint32_t id_block; // incremental id in the block container
	core::MetaEntry* meta;
	yon_gt* gt;
	yon_gt_summary* gt_sum;
	yon_occ* occ;
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
