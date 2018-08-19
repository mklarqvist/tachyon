#ifndef CORE_VARIANT_IMPORTER_STATS_H_
#define CORE_VARIANT_IMPORTER_STATS_H_

#include <cstdint>

namespace tachyon{
namespace support{

struct VariantImporterStats{
	typedef VariantImporterStats self_type;

	VariantImporterStats() :
		total_header_cost(0),
		total_gt_cost(0),
		total_ppa_cost(0),
		total_meta_cost(0),
		total_special_cost(0),
		total_info_cost(0),
		total_format_cost(0)
	{}

	inline const uint64_t getTotalCost(void) const{
		return(this->total_header_cost +
				this->total_gt_cost +
				this->total_ppa_cost +
				this->total_meta_cost +
				this->total_special_cost +
				this->total_info_cost +
				this->total_format_cost
			   );
	}


	uint64_t total_header_cost;
	uint64_t total_gt_cost;
	uint64_t total_ppa_cost;
	uint64_t total_meta_cost;
	uint64_t total_special_cost;
	uint64_t total_info_cost;
	uint64_t total_format_cost;
};

}
}



#endif /* CORE_VARIANT_IMPORTER_STATS_H_ */
