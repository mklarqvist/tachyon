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

bool VariantBlockContainer::readBlock(std::ifstream& stream, block_settings_type& settings){
	// Get info and format keys
	std::vector<U32> info_keys, format_keys;
	if(settings.info_ID_list.size()) info_keys = this->block_.intersectInfoKeys(settings.info_ID_list);
	else info_keys = this->block_.getInfoKeys();
	if(settings.format_ID_list.size()) format_keys = this->block_.intersectFormatKeys(settings.format_ID_list);
	else format_keys = this->block_.getFormatKeys();

	if(this->buildMapper(info_keys, format_keys, settings) == false)
		return false;

	// Todo: ascertain random access order is guaranteed

	if(settings.ppa.load){
		if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
			this->block_.ppa_manager.header = this->block_.footer.offset_ppa;
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offset_ppa.data_header.offset);
			stream >> this->block_.ppa_manager;
		}
	}

	if(settings.contig.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_contig, this->block_.meta_contig_container);
	}

	if(settings.positions.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_position, this->block_.meta_positions_container);
	}

	if(settings.controller.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_controllers, this->block_.meta_controller_container);
	}

	if(settings.quality.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_quality, this->block_.meta_quality_container);
	}

	if(settings.names.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_names, this->block_.meta_names_container);
	}

	if(settings.alleles.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_refalt, this->block_.meta_refalt_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_meta_alleles, this->block_.meta_alleles_container);
	}

	if(settings.genotypes_rle.load || settings.genotypes_all.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_8b, this->block_.gt_rle8_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_16b, this->block_.gt_rle16_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_32b, this->block_.gt_rle32_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_64b, this->block_.gt_rle64_container);
	}

	if(settings.genotypes_simple.load || settings.genotypes_all.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_simple8, this->block_.gt_simple8_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple16, this->block_.gt_simple16_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple32, this->block_.gt_simple32_container);
		this->block_.__loadContainer(stream, this->block_.footer.offset_gt_simple64, this->block_.gt_simple64_container);
	}

	if(settings.genotypes_support.load || settings.genotypes_all.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_gt_helper, this->block_.gt_support_data_container);
	}

	if(settings.set_membership.load || settings.genotypes_all.load){
		this->block_.__loadContainerSeek(stream, this->block_.footer.offset_meta_info_id, this->block_.meta_info_map_ids);
		this->block_.__loadContainer(stream, this->block_.footer.offset_meta_filter_id, this->block_.meta_filter_map_ids);
		this->block_.__loadContainer(stream, this->block_.footer.offset_meta_format_id, this->block_.meta_format_map_ids);
	}

	// Load all info
	if(settings.info_all.load && this->block_.footer.n_info_streams){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.__loadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
			this->mapper_.info_container_loaded_.at(i)(i, i, this->block_.footer.info_offsets[i].data_header.global_key, &this->block_.footer.info_offsets[i]);
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.info_ID_list.size()){
		this->mapper_.info_container_loaded_.resize(info_keys.size());
		// Ascertain that random access is linearly forward
		for(U32 i = 0; i < info_keys.size(); ++i){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local].data_header.offset);
			if(!stream.good()){
				std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for INFO field in block!" << std::endl;
				return false;
			}

			this->block_.__loadContainer(stream, this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local], this->block_.info_containers[i]);
			this->mapper_.info_container_loaded_.at(i)(i, this->mapper_.info_container_global_[info_keys[i]].stream_id_local, info_keys[i], &this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local]);
		}
	} // end case load_info_ID

	// Load all FORMAT data
	if(settings.format_all.load && this->block_.footer.n_format_streams){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.__loadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
			this->mapper_.format_container_loaded_.at(i)(i, i, this->block_.footer.format_offsets[i].data_header.global_key, &this->block_.footer.format_offsets[i]);
		}
		assert(this->block_.end_compressed_data_ == (U64)stream.tellg());
	} // If we have supplied a list of identifiers
	else if(settings.format_ID_list.size()){
		this->mapper_.format_container_loaded_.resize(format_keys.size());
		// Ascertain that random access is linearly forward
		for(U32 i = 0; i < format_keys.size(); ++i){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local].data_header.offset);
			if(!stream.good()){
				std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for FORMAT field in block!" << std::endl;
				return false;
			}

			this->block_.__loadContainer(stream, this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local], this->block_.format_containers[i]);
			this->mapper_.format_container_loaded_.at(i)(i, this->mapper_.format_container_global_[format_keys[i]].stream_id_local, format_keys[i], &this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local]);
		}
	} // end case load_info_ID

	stream.seekg(this->block_.end_block_); // seek to end-of-block
	return(true);
}

}
}
