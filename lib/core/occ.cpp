#include "occ.h"

namespace tachyon{

bool yon_occ::ReadTable(const std::string file_name,
                        const VariantHeader& header,
                        const char delimiter)
{
	std::ifstream f;
	f.open(file_name);
	if(!f.good()){
		std::cerr << "stream is bad: cannot open " << file_name << std::endl;
		return false;
	}

	uint32_t n_line = 0;
	std::string line;
	// Iterate over available lines in the input file.
	while(getline(f, line)){
		//std::cerr << line << std::endl;
		if(n_line++ == 0) continue;

		// Tokenize string with delimiter.
		std::vector<std::string> params = utility::split(line, delimiter, false);

		// Assert that the first column is an existing sample name.
		const int32_t sample_id = header.GetSampleId(params[0]);
		if(sample_id < 0){
			//std::cerr << "cannot find sample: " << params[0] << std::endl;
			//return false;
			continue;
		}

		// Iterate over tokens.
		for(U32 i = 1; i < params.size(); ++i){
			map_type::const_iterator it = this->map.find(params[i]);
			if(it == this->map.end()){
				// Not already set
				this->table.push_back(std::vector<uint32_t>(header.GetNumberSamples() + 1, 0));
				this->table[this->row_names.size()][sample_id + 1] = true;
				this->map[params[i]] = this->row_names.size();
				this->row_names.push_back(params[i]);
				std::cerr << "Adding group: " << params[i] << std::endl;
			} else {
				// Already set
				this->table[it->second][sample_id + 1] = 1;
			}
		}
	}

	/*
	for(U32 i = 0; i < this->table.size(); ++i){
		std::cerr << "Row " << i << "/" << this->table.size() << ": " << this->row_names[i];
		for(U32 j = 0; j < this->table[i].size(); ++j){
			std::cerr << "," << this->table[i][j];
		}
		std::cerr << std::endl;
	}
	*/

	return true;
}

bool yon_occ::BuildTable(void){
	this->occ.clear();
	if(this->table.size() == 0){
		return false;
	}

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	for(U32 i = 0; i < this->table.size(); ++i){
		//std::cerr << "Row " << i << "/" << this->table.size() << ": " << this->row_names[i];
		for(U32 j = 1; j < this->occ[i].size(); ++j){
			//std::cerr << "," << this->table[i][j];
			this->occ[i][j] += this->occ[i][j-1] + this->table[i][j];
			//std::cerr << "," << this->occ[i][j];
		}
		//std::cerr << std::endl;
	}
	return true;
}

bool yon_occ::BuildTable(const yon_gt_ppa& ppa){
	this->occ.clear();
	if(this->table.size() == 0){
		return false;
	}

	assert(ppa.n_samples + 1 == this->table[0].size());

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	for(U32 i = 0; i < this->table.size(); ++i){
		//std::cerr << "Row " << i << "/" << this->table.size() << ": " << this->row_names[i];
		for(U32 j = 1; j < this->occ[i].size(); ++j){
			//std::cerr << "," << this->table[i][j];
			this->occ[i][j] += this->occ[i][j - 1] + this->table[i][ppa[j - 1]];
			//std::cerr << "," << this->occ[i][j];
		}
		//std::cerr << std::endl;
		assert(this->occ[i][0] == 0);
	}
	return true;
}

}
