#include "variant_block_container.h"

namespace tachyon{
namespace containers{

const std::vector<bool> VariantBlockContainer::get_info_field_pattern_matches(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	int info_field_global_id = -2;
	if(this->header_->getInfoField(field_name, match)){
		info_field_global_id = match->IDX;
	}

	std::vector<bool> ret;
	if(info_field_global_id >= 0){
		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		S32 local_key = -1;
		for(U32 i = 0; i < this->getBlock().footer.n_info_streams; ++i){
			if(this->getBlock().footer.info_offsets[i].data_header.global_key == info_field_global_id){
				local_key = i;
			}
		}

		if(local_key == -1){
			//std::cerr << "could not find local" << std::endl;
			return ret;
		}

		ret.resize(this->getBlock().footer.n_info_patterns, false);
		for(U32 i = 0; i < this->getBlock().footer.n_info_patterns; ++i){
			//std::cerr << i << '\t' << this->getBlock().footer.info_bit_vectors[i][local_key] << std::endl;
			ret[i] = this->getBlock().footer.info_bit_vectors[i][local_key];
		}
	}
	return(ret);
}

const std::vector<bool> VariantBlockContainer::get_format_field_pattern_matches(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	int format_field_global_id = -2;
	if(this->header_->getFormatField(field_name, match)){
		format_field_global_id = match->IDX;
	}

	std::vector<bool> ret;
	if(format_field_global_id >= 0){
		S32 local_key = -1;
		for(U32 i = 0; i < this->getBlock().footer.n_format_streams; ++i){
			if(this->getBlock().footer.format_offsets[i].data_header.global_key == format_field_global_id){
				local_key = i;
			}
		}

		if(local_key == -1){
			//std::cerr << "could not find local" << std::endl;
			return ret;
		}

		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		ret.resize(this->getBlock().footer.n_format_patterns, false);
		for(U32 i = 0; i < this->getBlock().footer.n_format_patterns; ++i){
			//std::cerr << i << '\t' << this->getBlock().index_entry.format_bit_vectors[i][local_format_field_id] << std::endl;
			ret[i] = this->getBlock().footer.format_bit_vectors[i][local_key];
		}
	}
	return(ret);
}

}
}
