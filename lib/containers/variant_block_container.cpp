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
			this->block_.ppa_manager.header = this->block_.footer.offsets[YON_BLK_PPA];
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offsets[YON_BLK_PPA].data_header.offset);
			stream >> this->block_.ppa_manager;
		}
	}

	if(settings.contig.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_CONTIG], this->block_.base_containers[YON_BLK_CONTIG]);
	}

	if(settings.positions.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_POSITION], this->block_.base_containers[YON_BLK_POSITION]);
	}

	if(settings.controller.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_CONTROLLER], this->block_.base_containers[YON_BLK_CONTROLLER]);
	}

	if(settings.quality.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_QUALITY], this->block_.base_containers[YON_BLK_QUALITY]);
	}

	if(settings.names.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_NAMES], this->block_.base_containers[YON_BLK_NAMES]);
	}

	if(settings.alleles.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_REFALT], this->block_.base_containers[YON_BLK_REFALT]);
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_ALLELES], this->block_.base_containers[YON_BLK_ALLELES]);
	}

	if(settings.genotypes_rle.load || settings.genotypes_all.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_INT8], this->block_.base_containers[YON_BLK_GT_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT16], this->block_.base_containers[YON_BLK_GT_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT32], this->block_.base_containers[YON_BLK_GT_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT64], this->block_.base_containers[YON_BLK_GT_INT64]);
	}

	if(settings.genotypes_simple.load || settings.genotypes_all.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT8], this->block_.base_containers[YON_BLK_GT_S_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT16], this->block_.base_containers[YON_BLK_GT_S_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT32], this->block_.base_containers[YON_BLK_GT_S_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT64], this->block_.base_containers[YON_BLK_GT_S_INT64]);
	}

	if(settings.genotypes_support.load || settings.genotypes_all.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_SUPPORT], this->block_.base_containers[YON_BLK_GT_SUPPORT]);
	}

	if(settings.set_membership.load || settings.genotypes_all.load){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_ID_INFO], this->block_.base_containers[YON_BLK_ID_INFO]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_ID_FORMAT], this->block_.base_containers[YON_BLK_ID_FORMAT]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_ID_FILTER], this->block_.base_containers[YON_BLK_ID_FILTER]);
	}

	// Load all info
	if(settings.info_all.load && this->block_.footer.n_info_streams){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
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

			this->block_.LoadContainer(stream, this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local], this->block_.info_containers[i]);
			this->mapper_.info_container_loaded_.at(i)(i, this->mapper_.info_container_global_[info_keys[i]].stream_id_local, info_keys[i], &this->block_.footer.info_offsets[this->mapper_.info_container_global_[info_keys[i]].stream_id_local]);
		}
	} // end case load_info_ID

	// Load all FORMAT data
	if(settings.format_all.load && this->block_.footer.n_format_streams){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
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

			this->block_.LoadContainer(stream, this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local], this->block_.format_containers[i]);
			this->mapper_.format_container_loaded_.at(i)(i, this->mapper_.format_container_global_[format_keys[i]].stream_id_local, format_keys[i], &this->block_.footer.format_offsets[this->mapper_.format_container_global_[format_keys[i]].stream_id_local]);
		}
	} // end case load_info_ID

	stream.seekg(this->block_.end_block_); // seek to end-of-block
	return(true);
}

VariantReaderObjects& VariantBlockContainer::loadObjects(objects_type& objects, block_settings_type& block_settings) const{
	// New meta container
	objects.meta_container = new meta_container_type(this->getBlock());

	// New genotype containers if aplicable
	if(this->getBlock().header.controller.hasGT && block_settings.genotypes_all.load){
		objects.genotype_container = new gt_container_type(this->getBlock(), *objects.meta_container);
		objects.genotype_summary   = new objects_type::genotype_summary_type(10);
	}

	// FORMAT-specific containers
	// Store as double pointers to avoid memory collisions because
	// FORMAT containers have different intrinsic class members
	objects.n_loaded_format   = this->getMapper().getNumberFormatLoaded();
	objects.format_containers = new format_interface_type*[objects.n_loaded_format];

	if(objects.n_loaded_format){
		for(U32 i = 0; i < objects.n_loaded_format; ++i){
			const U32 global_key = this->getMapper().getLoadedFormat(i).stream_id_global;

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->get_format_field_pattern_matches(this->header_->format_fields[global_key].ID);

			if(this->header_->format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.format_containers[i] = new containers::FormatContainer<S32>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->getSampleNumber());
				objects.format_field_names.push_back(this->header_->format_fields[global_key].ID);
			} else if(this->header_->format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header_->format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.format_containers[i] = new containers::FormatContainer<std::string>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->getSampleNumber());
				objects.format_field_names.push_back(this->header_->format_fields[global_key].ID);
			} else if(this->header_->format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.format_containers[i] = new containers::FormatContainer<float>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->getSampleNumber());
				objects.format_field_names.push_back(this->header_->format_fields[global_key].ID);
			} else {
				objects.format_containers[i] = new containers::FormatContainer<U32>;
				objects.format_field_names.push_back(this->header_->format_fields[global_key].ID);
			}
			objects.format_container_map[this->header_->format_fields[global_key].ID] = objects.format_containers[i];
		}
	}

	// INFO-specific containers
	// Store as double pointers to avoid memory collisions because
	// INFO containers have different class members
	objects.n_loaded_info   = this->getMapper().getNumberInfoLoaded();
	objects.info_containers = new info_interface_type*[objects.n_loaded_info];

	if(objects.n_loaded_info){
		for(U32 i = 0; i < objects.n_loaded_info; ++i){
			const U32 global_key = this->getMapper().getLoadedInfo(i).stream_id_global;

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->get_info_field_pattern_matches(this->header_->info_fields[global_key].ID);

			if(this->header_->info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.info_containers[i] = new containers::InfoContainer<S32>(this->getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header_->info_fields[global_key].ID);
			} else if(this->header_->info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header_->info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.info_containers[i] = new containers::InfoContainer<std::string>(this->getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header_->info_fields[global_key].ID);
			} else if(this->header_->info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.info_containers[i] = new containers::InfoContainer<float>(this->getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header_->info_fields[global_key].ID);
			} else {
				objects.info_containers[i] = new containers::InfoContainer<U32>();
				objects.info_field_names.push_back(this->header_->info_fields[global_key].ID);
			}
			objects.info_container_map[this->header_->info_fields[global_key].ID] = objects.info_containers[i];
		}
	}


	// If we want to drop records that do not have all/any of the fields we desire
	// then we create a vector of size N_PATTERNS and set those that MATCH to TRUE
	// this allows for filtering in O(1)-time
	//
	// This vector stores the number of INFO fields having set membership with this
	// particular hash pattern
	objects.info_id_fields_keep   = std::vector<U32>(this->getBlock().footer.n_info_patterns,   0);
	objects.format_id_fields_keep = std::vector<U32>(this->getBlock().footer.n_format_patterns, 0);

	// This vector of vectors keeps the local INFO identifiers for the matched global
	// identifiers in a given hash pattern
	//
	// For example: x[5] contains the local IDs for loaded INFO streams for pattern ID 5
	objects.local_match_keychain_info   = std::vector< std::vector<U32> >(this->getBlock().footer.n_info_patterns);
	objects.local_match_keychain_format = std::vector< std::vector<U32> >(this->getBlock().footer.n_format_patterns);

	// If loading all INFO values then return them in the ORIGINAL order
	if(block_settings.info_all.load){
		for(U32 i = 0; i < this->getBlock().footer.n_info_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->getBlock().footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < objects.n_loaded_info; ++k){ // Number of loaded INFO identifiers
					// Global
					if(this->getBlock().footer.info_offsets[this->getBlock().footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == this->getMapper().getLoadedInfo(k).offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_id_fields_keep[i];
					}
				}
			}
		}
	}
	// If loading custom INFO fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->getBlock().footer.n_info_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < objects.n_loaded_info; ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->getBlock().footer.info_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->getBlock().footer.info_offsets[this->getBlock().footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == this->getMapper().getLoadedInfo(k).offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_id_fields_keep[i];
					}
				}
			}
		}
	}

	// For FORMAT
	// If loading all FORMAT values then return them in the ORIGINAL order
	if(block_settings.format_all.load){
		for(U32 i = 0; i < this->getBlock().footer.n_format_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->getBlock().footer.format_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < objects.n_loaded_format; ++k){ // Number of loaded INFO identifiers
					if(this->getBlock().footer.format_offsets[this->getBlock().footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == this->getMapper().getLoadedFormat(k).offset->data_header.global_key){
						objects.local_match_keychain_format[i].push_back(k);
						++objects.format_id_fields_keep[i];
					}
				}
			}
		}
	}
	// If loading custom FORMAT fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->getBlock().footer.n_format_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < objects.n_loaded_format; ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->getBlock().footer.format_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->getBlock().footer.format_offsets[this->getBlock().footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == this->getMapper().getLoadedFormat(k).offset->data_header.global_key){
						objects.local_match_keychain_format[i].push_back(k);
						++objects.format_id_fields_keep[i];
					}
				}
			}
		}
	}


	// If we want to compute additional genotypic summary statistics (triggered by -X flag)
	// then we need to make sure we don't accidentally add annotations to fields that already
	// exists (e.g. if INFO=AC already exists)
	//
	// Preprocessing step:
	// Cycle over INFO patterns and see if any of the custom FIELDS are set
	// FS_A, AN, NM, NPM, AC, AC_FW, AC_REV, AF, HWE_P, VT, MULTI_ALLELIC
	std::vector<std::string> ADDITIONAL_INFO = {"FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"};
	U16 execute_mask = 0;

	// Step 1: Find INFO
	std::vector< std::pair<U32, U32> > additional_local_keys_found;
	for(U32 i = 0; i < ADDITIONAL_INFO.size(); ++i){
		if(this->header_->has_info_field(ADDITIONAL_INFO[i])){
			const core::HeaderMapEntry* map = this->header_->getInfoField(ADDITIONAL_INFO[i]);
			// Find local key
			for(U32 k = 0; k < objects.n_loaded_info; ++k){
				if(this->getBlock().info_containers[k].header.getGlobalKey() == map->IDX){
					execute_mask |= 1 << i;
					additional_local_keys_found.push_back(std::pair<U32,U32>(k,i));
				}
			}
		}
	}

	// Step 2: Cycle over patterns to find existing INFO fields
	// Cycle over INFO patterns
	objects.additional_info_execute_flag_set = std::vector< U16 >(1, 65535);
	if(ADDITIONAL_INFO.size()){
		objects.additional_info_execute_flag_set.reserve(this->getBlock().footer.n_info_patterns);
		for(U32 i = 0; i < this->getBlock().footer.n_info_patterns; ++i){
			objects.additional_info_execute_flag_set[i] = (1 << ADDITIONAL_INFO.size()) - 1;
			for(U32 j = 0; j < additional_local_keys_found.size(); ++j){
				if(this->getBlock().footer.info_bit_vectors[i][j]){
					objects.additional_info_execute_flag_set[i] &= ~(1 << additional_local_keys_found[j].second);
				}
			}
		}
	}

	//

	return(objects);
}

}
}
