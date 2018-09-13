#include "occ.h"

namespace tachyon{

bool yon_occ::ReadTable(const std::string file_name,
                        const VariantHeader& header,
                        const char delimiter)
{
	std::ifstream f;
	f.open(file_name);
	if(!f.good()){
		std::cerr << utility::timestamp("ERROR") << "Stream is bad! Cannot open " << file_name << "..." << std::endl;
		return false;
	}

	uint32_t n_line = 0;
	std::string line;
	// Iterate over available lines in the input file.
	while(getline(f, line)){
		//std::cerr << line << std::endl;

		// Tokenize string with delimiter.
		std::vector<std::string> params = utility::split(line, delimiter, false);

		// Assert that the first column is an existing sample name.
		const int32_t sample_id = header.GetSampleId(params[0]);
		if(sample_id < 0){
			std::cerr << utility::timestamp("WARNING") << "Cannot find sample \"" << params[0] << "\" in groupings file..." << std::endl;
			continue;
		}

		// Iterate over tokens.
		for(uint32_t i = 1; i < params.size(); ++i){
			map_type::const_iterator it = this->map.find(params[i]);
			if(it == this->map.end()){
				// Not already set
				this->table.push_back(std::vector<uint32_t>(header.GetNumberSamples() + 1, 0));
				this->table.back()[sample_id + 1] = true;
				this->map[params[i]] = this->row_names.size();
				this->row_names.push_back(params[i]);
				//std::cerr << "Adding group: " << params[i] << " for " << this->table.back().size() << " samples" << std::endl;
			} else {
				// Already set
				this->table[it->second][sample_id + 1] = true;
			}
		}
	}

	return true;
}

bool yon_occ::BuildTable(void){
	this->occ.clear();
	if(this->table.size() == 0){
		return false;
	}

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	this->cum_sums = std::vector< uint32_t >( this->occ.size() );

	for(uint32_t i = 0; i < this->table.size(); ++i){
		assert(this->table[i][0] == 0);
		for(uint32_t j = 1; j < this->occ[i].size(); ++j)
			this->occ[i][j] += this->occ[i][j-1] + this->table[i][j];

		this->cum_sums[i] = this->occ[i].back();
	}

	return true;
}

bool yon_occ::BuildTable(const yon_gt_ppa& ppa){
	this->occ.clear();
	if(this->table.size() == 0){
		return false;
	}

	assert(ppa.n_s + 1 == this->table[0].size());

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	this->cum_sums = std::vector< uint32_t >( this->occ.size() );

	for(uint32_t i = 0; i < this->table.size(); ++i){
		assert(this->table[i][0] == 0);
		for(uint32_t j = 1; j < this->occ[i].size(); ++j)
			this->occ[i][j] += this->occ[i][j - 1] + this->table[i][ppa[j - 1] + 1];

		assert(this->occ[i][0] == 0);
		this->cum_sums[i] = this->occ[i].back();
	}

	return true;
}

}
