#ifndef CORE_IMPORTERSTATS_H_
#define CORE_IMPORTERSTATS_H_

namespace Tachyon{
namespace Support{

struct ImporterStats{
	typedef ImporterStats self_type;

	ImporterStats() :
		total_header_cost(0),
		total_gt_cost(0),
		total_ppa_cost(0),
		total_meta_cost(0),
		total_info_cost(0),
		total_format_cost(0),
		total_rest_cost(0)
	{}

	inline const U64 getTotalCost(void) const{
		return(this->total_header_cost +
				this->total_gt_cost +
				this->total_ppa_cost +
				this->total_meta_cost +
				this->total_info_cost +
				this->total_format_cost +
				this->total_rest_cost
			   );
	}

	// Used for debugging only
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << "Header: " << entry.total_header_cost << '\n'
			<< "GT: " << entry.total_gt_cost << '\n'
			<< "PPA: " << entry.total_ppa_cost << '\n'
			<< "Meta: " << entry.total_meta_cost << '\n'
			<< "INFO: " << entry.total_info_cost << '\n'
			<< "FORMAT: " << entry.total_format_cost << '\n'
			<< "Residual: " << entry.total_rest_cost;

		return(out);
	}


	U64 total_header_cost;
	U64 total_gt_cost;
	U64 total_ppa_cost;
	U64 total_meta_cost;
	U64 total_info_cost;
	U64 total_format_cost;
	U64 total_rest_cost;
};

}
}



#endif /* CORE_IMPORTERSTATS_H_ */
