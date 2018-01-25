#ifndef CORE_IMPORTERSTATS_H_
#define CORE_IMPORTERSTATS_H_

namespace tachyon{
namespace support{

struct ImporterStats{
	typedef ImporterStats self_type;

	ImporterStats() :
		total_header_cost(0),
		total_gt_cost(0),
		total_ppa_cost(0),
		total_meta_cost(0),
		total_special_cost(0),
		total_info_cost(0),
		total_format_cost(0)
	{}

	inline const U64 getTotalCost(void) const{
		return(this->total_header_cost +
				this->total_gt_cost +
				this->total_ppa_cost +
				this->total_meta_cost +
				this->total_special_cost +
				this->total_info_cost +
				this->total_format_cost
			   );
	}


	U64 total_header_cost;
	U64 total_gt_cost;
	U64 total_ppa_cost;
	U64 total_meta_cost;
	U64 total_special_cost;
	U64 total_info_cost;
	U64 total_format_cost;
};

}
}



#endif /* CORE_IMPORTERSTATS_H_ */
