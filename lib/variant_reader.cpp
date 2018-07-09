#include "variant_reader.h"

namespace tachyon{

VariantReader::VariantReader() :
	filesize(0)
{}

VariantReader::VariantReader(const std::string& filename) :
	filesize(0)
{}

VariantReader::~VariantReader(){}

VariantReader::VariantReader(const self_type& other) :
	filesize(other.filesize),
	block_settings(other.block_settings),
	settings(other.settings),
	header(other.header),
	footer(other.footer),
	index(other.index),
	checksums(other.checksums),
	keychain(other.keychain)
{
	this->stream.open(this->settings.input, std::ios::in | std::ios::binary);
}

bool VariantReader::open(void){
	if(this->settings.input.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	this->stream.open(this->settings.input, std::ios::binary | std::ios::in | std::ios::ate);
	this->filesize = (U64)this->stream.tellg();

	if(this->filesize <= YON_FOOTER_LENGTH){
		std::cerr << utility::timestamp("ERROR") << "File is corrupted!" << std::endl;
		return false;
	}

	// Seek to start of footer
	this->stream.seekg(this->filesize - YON_FOOTER_LENGTH);
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to seek in file!" << std::endl;
		return false;
	}
	this->stream >> this->footer;

	// Validate footer
	if(this->footer.validate() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate footer!" << std::endl;
		return false;
	}

	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to read file!" << std::endl;
		return false;
	}

	// Seek to start of file
	this->stream.seekg(0);
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to rewind file!" << std::endl;
		return false;
	}

	// Load header
	//this->stream >> this->header;
	char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	this->stream.read(&magic_string[0], tachyon::constants::FILE_HEADER_LENGTH);
	if(strncmp(&magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate Tachyon magic string!" << std::endl;
		return false;
	}

	containers::DataContainer header_container;
	this->stream >> header_container.header;
	header_container.buffer_data.resize(header_container.header.data_header.cLength);
	this->stream.read(header_container.buffer_data.data(), header_container.header.data_header.cLength);
	header_container.buffer_data.n_chars = header_container.header.data_header.cLength;
	if(!this->codec_manager.zstd_codec.decompress(header_container)){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	header_container.buffer_data_uncompressed >> this->header; // parse header from buffer

	if(!this->header.header_magic.validate()){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate header!" << std::endl;
		return false;
	}

	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to get header!" << std::endl;
		return false;
	}

	this->variant_container << this->header;

	// Keep track of start position
	const U64 return_pos = this->stream.tellg();
	this->stream.seekg(this->footer.offset_end_of_data);
	this->stream >> this->index;
	this->stream >> this->checksums;
	this->stream.seekg(return_pos);
	return(this->stream.good());
}

bool VariantReader::nextBlock(){
	// If the stream is faulty then return
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->stream.tellg() == this->footer.offset_end_of_data)
		return false;

	// Reset and re-use
	this->variant_container.reset();

	if(!this->variant_container.getBlock().readHeaderFooter(this->stream))
		return false;

	if(!this->codec_manager.zstd_codec.decompress(this->variant_container.getBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	this->variant_container.getBlock().footer_support.buffer_data_uncompressed >> this->variant_container.getBlock().footer;

	// Parse settings
	this->parseSettings();
	if(this->variant_container.build() == false){
		return false;
	}

	// Attempts to read a YON block with the provided
	if(!this->variant_container.getBlock().read(this->stream, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.anyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption::EncryptionDecorator e;
		if(!e.decryptAES256(this->variant_container.getBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->variant_container.getBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

bool VariantReader::getBlock(const index_entry_type& index_entry){
	// If the stream is faulty then return
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// Reset and re-use
	this->variant_container.getBlock().clear();

	// Seek to target block ID with the help of the index
	this->stream.seekg(index_entry.byte_offset);
	if(this->stream.good() == false){
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to seek to given offset using target index entry!" << std::endl;
		return(false);
	}

	if(!this->variant_container.getBlock().readHeaderFooter(this->stream))
		return false;

	if(!this->codec_manager.zstd_codec.decompress(this->variant_container.getBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
		return(false);
	}
	this->variant_container.getBlock().footer_support.buffer_data_uncompressed >> this->variant_container.getBlock().footer;
	this->parseSettings();

	// Attempts to read a YON block with the provided
	if(!this->variant_container.getBlock().read(this->stream, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.getBlock().header.controller.anyEncrypted){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption::EncryptionDecorator e;
		if(!e.decryptAES256(this->variant_container.getBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->variant_container.getBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

VariantReader::block_entry_type VariantReader::getBlock(){
	// If the stream is faulty then return
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR") << "Corrupted!" << std::endl;
		return block_entry_type();
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->stream.tellg() == this->footer.offset_end_of_data)
		return block_entry_type();


	block_entry_type block;
	const size_t position = this->stream.tellg();

	// Attempts to read a YON block with the provided
	// settings
	if(!block.read(this->stream, this->block_settings)){
		this->stream.seekg(position);
		return block;
	}

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->variant_container.getBlock())){
		this->stream.seekg(position);
		return block;
	}

	this->stream.seekg(position);

	// All passed
	return block;
}

const int VariantReader::has_format_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getFormatField(field_name, match)){
		return(match->IDX);
	}
	return(-2);
}

const int VariantReader::has_info_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getInfoField(field_name, match)){
		return(match->IDX);
	}
	return(-2);
}

const int VariantReader::has_filter_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getFilterField(field_name, match)){
		return(match->IDX);
	}
	return(-2);
}

VariantReaderObjects& VariantReader::loadObjects(objects_type& objects) const{
	// New meta container
	objects.meta_container = new meta_container_type(this->variant_container.getBlock());

	// New genotype containers if aplicable
	if(this->variant_container.getBlock().header.controller.hasGT && block_settings.genotypes_all.load){
		objects.genotype_container = new gt_container_type(this->variant_container.getBlock(), *objects.meta_container);
		objects.genotype_summary   = new objects_type::genotype_summary_type(10);
	}

	// FORMAT-specific containers
	// Store as double pointers to avoid memory collisions because
	// FORMAT containers have different class members
	objects.n_loaded_format   = this->variant_container.getBlock().n_format_loaded;
	objects.format_containers = new format_interface_type*[objects.n_loaded_format];

	if(this->variant_container.getBlock().n_format_loaded){
		for(U32 i = 0; i < this->variant_container.getBlock().n_format_loaded; ++i){
			// Remap GLOBAL identifier to the LOCAL identifier
			const U32 global_key = block_settings.load_format_ID_loaded[i].offset->data_header.global_key;
			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->get_format_field_pattern_matches(this->header.format_fields[global_key].ID);

			if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.format_containers[i] = new containers::FormatContainer<S32>(this->variant_container.getBlock().format_containers[i], *objects.meta_container, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header.format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.format_containers[i] = new containers::FormatContainer<std::string>(this->variant_container.getBlock().format_containers[i], *objects.meta_container, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.format_containers[i] = new containers::FormatContainer<float>(this->variant_container.getBlock().format_containers[i], *objects.meta_container, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else {
				objects.format_containers[i] = new containers::FormatContainer<U32>;
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			}
		}
	}

	// INFO-specific containers
	// Store as double pointers to avoid memory collisions because
	// INFO containers have different class members
	objects.n_loaded_info   = this->variant_container.getBlock().n_info_loaded;
	objects.info_containers = new info_interface_type*[objects.n_loaded_info];

	if(this->variant_container.getBlock().n_info_loaded){
		for(U32 i = 0; i < this->variant_container.getBlock().n_info_loaded; ++i){
			// Remap GLOBAL identifier to the LOCAL identifier
			const U32 global_key = block_settings.load_info_ID_loaded[i].offset->data_header.global_key;
			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[global_key].ID);

			if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.info_containers[i] = new containers::InfoContainer<S32>(this->variant_container.getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header.info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.info_containers[i] = new containers::InfoContainer<std::string>(this->variant_container.getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.info_containers[i] = new containers::InfoContainer<float>(this->variant_container.getBlock().info_containers[i], *objects.meta_container, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else {
				objects.info_containers[i] = new containers::InfoContainer<U32>();
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			}
		}
	}

	// If we want to drop records that do not have all/any of the fields we desire
	// then we create a vector of size N_PATTERNS and set those that MATCH to TRUE
	// this allows for filtering in O(1)-time
	//
	// This vector stores the number of INFO fields having set membership with this
	// particular hash pattern
	objects.info_id_fields_keep   = std::vector<U32>(this->variant_container.getBlock().footer.n_info_patterns,   0);
	objects.format_id_fields_keep = std::vector<U32>(this->variant_container.getBlock().footer.n_format_patterns, 0);

	// This vector of vectors keeps the local INFO identifiers for the matched global
	// identifiers in a given hash pattern
	//
	// For example: x[5] contains the local IDs for loaded INFO streams for pattern ID 5
	objects.local_match_keychain_info   = std::vector< std::vector<U32> >(this->variant_container.getBlock().footer.n_info_patterns);
	objects.local_match_keychain_format = std::vector< std::vector<U32> >(this->variant_container.getBlock().footer.n_format_patterns);

	// If loading all INFO values then return them in the ORIGINAL order
	if(this->block_settings.info_all.load){
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->variant_container.getBlock().footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < block_settings.load_info_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
					if(this->variant_container.getBlock().footer.info_offsets[this->variant_container.getBlock().footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == block_settings.load_info_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_id_fields_keep[i];
					}
				}
			}
		}
	}
	// If loading custom INFO fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < block_settings.load_info_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->variant_container.getBlock().footer.info_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->variant_container.getBlock().footer.info_offsets[this->variant_container.getBlock().footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == block_settings.load_info_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_id_fields_keep[i];
					}
				}
			}
		}
	}

	// For FORMAT
	// If loading all FORMAT values then return them in the ORIGINAL order
	if(this->block_settings.format_all.load){
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->variant_container.getBlock().footer.format_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < block_settings.load_format_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
					if(this->variant_container.getBlock().footer.format_offsets[this->variant_container.getBlock().footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == block_settings.load_format_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_format[i].push_back(k);
						++objects.format_id_fields_keep[i];
					}
				}
			}
		}
	}
	// If loading custom FORMAT fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < block_settings.load_format_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->variant_container.getBlock().footer.format_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->variant_container.getBlock().footer.format_offsets[this->variant_container.getBlock().footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == block_settings.load_format_ID_loaded[k].offset->data_header.global_key){
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
		if(this->header.has_info_field(ADDITIONAL_INFO[i])){
			const core::HeaderMapEntry* map = this->header.getInfoField(ADDITIONAL_INFO[i]);
			// Find local key
			for(U32 k = 0; k < this->block_settings.load_info_ID_loaded.size(); ++k){
				if(this->variant_container.getBlock().info_containers[k].header.getGlobalKey() == map->IDX){
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
		objects.additional_info_execute_flag_set.reserve(this->variant_container.getBlock().footer.n_info_patterns);
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_patterns; ++i){
			objects.additional_info_execute_flag_set[i] = (1 << ADDITIONAL_INFO.size()) - 1;
			for(U32 j = 0; j < additional_local_keys_found.size(); ++j){
				if(this->variant_container.getBlock().footer.info_bit_vectors[i][j]){
					objects.additional_info_execute_flag_set[i] &= ~(1 << additional_local_keys_found[j].second);
				}
			}
		}
	}

	//

	return(objects);
}

const std::vector<bool> VariantReader::get_info_field_pattern_matches(const std::string& field_name) const{
	int global_info_value = this->has_info_field(field_name);

	std::vector<bool> ret;
	if(global_info_value >= 0){
		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		S32 local_key = -1;
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_streams; ++i){
			if(this->variant_container.getBlock().footer.info_offsets[i].data_header.global_key == global_info_value){
				local_key = i;
			}
		}

		if(local_key == -1){
			//std::cerr << "could not find local" << std::endl;
			return ret;
		}

		ret.resize(this->variant_container.getBlock().footer.n_info_patterns, false);
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_patterns; ++i){
			//std::cerr << i << '\t' << this->variant_container.getBlock().footer.info_bit_vectors[i][local_key] << std::endl;
			ret[i] = this->variant_container.getBlock().footer.info_bit_vectors[i][local_key];
		}
	}
	return(ret);
}

const std::vector<bool> VariantReader::get_format_field_pattern_matches(const std::string& field_name) const{
	int global_format_value = this->has_format_field(field_name);
	std::vector<bool> ret;
	if(global_format_value >= 0){
		S32 local_key = -1;
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_streams; ++i){
			if(this->variant_container.getBlock().footer.format_offsets[i].data_header.global_key == global_format_value){
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
		ret.resize(this->variant_container.getBlock().footer.n_format_patterns, false);
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_patterns; ++i){
			//std::cerr << i << '\t' << this->variant_container.getBlock().index_entry.format_bit_vectors[i][local_format_field_id] << std::endl;
			ret[i] = this->variant_container.getBlock().footer.format_bit_vectors[i][local_key];
		}
	}
	return(ret);
}

bool VariantReader::parseSettings(void){
	this->block_settings.load_info_ID_loaded.clear();
	this->block_settings.load_format_ID_loaded.clear();
	this->block_settings.format_map.clear();
	this->block_settings.info_map.clear();

	this->block_settings.format_map.resize(this->header.header_magic.n_format_values);
	this->block_settings.info_map.resize(this->header.header_magic.n_info_values);

	// Map INFO
	if(this->block_settings.info_all.load == false){ // prevent double load
		U32 info_matches = 0;
		for(U32 i = 0; i < this->block_settings.info_list.size(); ++i){
			const S32 global_key = this->has_info_field(this->block_settings.info_list[i]);
			if(global_key >= 0){
				S32 local_key = -1;
				for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_streams; ++i){
					if(this->variant_container.getBlock().footer.info_offsets[i].data_header.global_key == global_key){
						local_key = i;
						this->block_settings.info_ID_list.push_back(local_key);

						this->block_settings.info_map[global_key].stream_id_local  = local_key;
						this->block_settings.info_map[global_key].stream_id_global = global_key;
						this->block_settings.info_map[global_key].load_order_index = info_matches;

						this->block_settings.load_info_ID_loaded.push_back(
													SettingsMap(
															info_matches++, // iterator value
															local_key,      // local index id
															&this->variant_container.getBlock().footer.info_offsets[local_key]) // offset
													);
						break;
					}
				}

				//if(local_key == -1){
					//std::cerr << "could not find local" << std::endl;
				//}
			}
		}
	} else {
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_info_streams; ++i){
			const U32 global_key = this->variant_container.getBlock().footer.info_offsets[i].data_header.global_key;
			this->block_settings.info_map[i].stream_id_local  = i;
			this->block_settings.info_map[i].stream_id_global = global_key;
			this->block_settings.info_map[i].load_order_index = i;
		}
	}

	// Map FORMAT
	if(this->block_settings.format_all.load == false){ // prevent double load
		U32 format_matches = 0;
		for(U32 i = 0; i < this->block_settings.format_list.size(); ++i){
			const S32 global_key = this->has_format_field(this->block_settings.format_list[i]);
			if(global_key >= 0){
				S32 local_key = -1;
				for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_streams; ++i){
					if(this->variant_container.getBlock().footer.format_offsets[i].data_header.global_key == global_key){
						local_key = i;
						this->block_settings.format_ID_list.push_back(local_key);

						this->block_settings.format_map[global_key].stream_id_local  = local_key;
						this->block_settings.format_map[global_key].stream_id_global = global_key;
						this->block_settings.format_map[global_key].load_order_index = format_matches;

						this->block_settings.load_format_ID_loaded.push_back(
													SettingsMap(
															format_matches++, // iterator value
															local_key,        // local index id
															&this->variant_container.getBlock().footer.format_offsets[local_key]) // offset
													);
						break;
					}
				}

				//if(local_key == -1){
					//std::cerr << "could not find local" << std::endl;
				//}
			}
		}
	} else {
		for(U32 i = 0; i < this->variant_container.getBlock().footer.n_format_streams; ++i){
			const U32 global_key = this->variant_container.getBlock().footer.format_offsets[i].data_header.global_key;
			this->block_settings.format_map[i].stream_id_local  = i;
			this->block_settings.format_map[i].stream_id_global = global_key;
			this->block_settings.format_map[i].load_order_index = i;
		}
	}

	return(true);
}


void VariantReader::printFILTER(buffer_type& outputBuffer,
				 const U32& position,
				 const objects_type& objects) const
{
	if(outputBuffer.back() != '\t') outputBuffer += '\t';

	if(this->block_settings.display_filter && this->variant_container.getBlock().footer.n_filter_streams){
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			// Local key -> global key
			outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
			}
		} else
			outputBuffer += '.';
	} else {
		outputBuffer += '.';
	}
}

void VariantReader::printFILTERCustom(buffer_type& outputBuffer,
				 const U32& position,
				 const objects_type& objects) const
{
	if(this->block_settings.display_filter && this->variant_container.getBlock().footer.n_filter_streams){
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			if(outputBuffer.back() != this->block_settings.custom_delimiter_char) outputBuffer += this->block_settings.custom_delimiter_char;

			// Local key -> global key
			outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
			}
		} else
			outputBuffer += '.';
	} else {
		outputBuffer += '.';
	}
}


void VariantReader::printFILTERJSON(buffer_type& outputBuffer,
					 const U32& position,
					 const objects_type& objects) const
{
	if(this->variant_container.getBlock().footer.n_filter_streams){
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->variant_container.getBlock().footer.filter_bit_vectors[(*objects.meta_container)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Local key -> global key
			outputBuffer += "\"FILTER-";
			outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			outputBuffer += "\":true";
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ',';
				outputBuffer += "\"FILTER-";
				outputBuffer += this->header.filter_fields[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
				outputBuffer += "\":true";
			}
		}
	}
}

void VariantReader::printFORMATVCF(buffer_type& buffer,
	                                const char& delimiter,
	                                 const U32& position,
	                        const objects_type& objects,
	               std::vector<core::GTObject>& genotypes_unpermuted) const
{
	if(this->block_settings.format_all.display && this->variant_container.getBlock().n_format_loaded){
		if(this->variant_container.getBlock().n_format_loaded){
			const U32& n_format_keys = this->variant_container.getBlock().footer.format_bit_vectors[objects.meta_container->at(position).format_pattern_id].n_keys;
			const U32* format_keys   = this->variant_container.getBlock().footer.format_bit_vectors[objects.meta_container->at(position).format_pattern_id].local_keys;
			if(n_format_keys){
				if(buffer.back() != delimiter) buffer += delimiter;

				// Print key map
				buffer += this->header.format_fields[this->variant_container.getBlock().footer.format_offsets[format_keys[0]].data_header.global_key].ID;
				for(U32 i = 1; i < n_format_keys; ++i){
					buffer += ':';
					buffer += this->header.format_fields[this->variant_container.getBlock().footer.format_offsets[format_keys[i]].data_header.global_key].ID;
				}
				buffer += delimiter;

				// Todo: print if no GT data
				// Begin print FORMAT data for each sample
				if(this->variant_container.getBlock().header.controller.hasGT){
					if(this->block_settings.ppa.load && this->variant_container.getBlock().header.controller.hasGTPermuted){
						objects.genotype_container->at(position).getObjects(this->header.getSampleNumber(), genotypes_unpermuted, this->variant_container.getBlock().ppa_manager);
					} else {
						objects.genotype_container->at(position).getObjects(this->header.getSampleNumber(), genotypes_unpermuted);
					}

					buffer << genotypes_unpermuted[0];
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						objects.format_containers[format_keys[i]]->to_vcf_string(buffer, position, 0);
					}

					for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
						buffer += delimiter;
						buffer << genotypes_unpermuted[s];
						for(U32 i = 1; i < n_format_keys; ++i){
							buffer  += ':';
							objects.format_containers[format_keys[i]]->to_vcf_string(buffer, position, s);
						}
					}
				} else { // No genotype data available
					objects.format_containers[format_keys[0]]->to_vcf_string(buffer, position, 0);
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						objects.format_containers[format_keys[i]]->to_vcf_string(buffer, position, 0);
					}

					for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
						buffer += delimiter;
						objects.format_containers[format_keys[0]]->to_vcf_string(buffer, position, s);
						for(U32 i = 1; i < n_format_keys; ++i){
							buffer  += ':';
							objects.format_containers[format_keys[i]]->to_vcf_string(buffer, position, s);
						}
					}
				}

			} else { // have no keys
				if(buffer.back() != delimiter) buffer += delimiter;
				buffer += ".";
				//buffer += delimiter;
			}
		}
	}
}

void VariantReader::printFORMATVCF(buffer_type& buffer, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const
{
	return(this->printFORMATVCF(buffer, '\t', position, objects, genotypes_unpermuted));
}

void VariantReader::printFORMATCustom(buffer_type& outputBuffer,
					   const char& delimiter,
					   const U32& position,
					   const objects_type& objects,
					   std::vector<core::GTObject>& genotypes_unpermuted) const
{
	if(block_settings.format_all.display || this->variant_container.getBlock().n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta_container->at(position).format_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// Print key map
			outputBuffer += objects.format_field_names[targetKeys[0]];
			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ':';
				outputBuffer += objects.format_field_names[targetKeys[i]];
			};

			outputBuffer += delimiter;

			// First individual
			//if(this->header.getSampleNumber() > 1) outputBuffer += '[';

			//if(targetKeys.size() > 1) outputBuffer += '[';
			objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ':';
				objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
			}
			//if(targetKeys.size() > 1) outputBuffer += ']';

			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += delimiter;
				//if(targetKeys.size() > 1) outputBuffer += '[';
				objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer  += ':';
					objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
				}
				//if(targetKeys.size() > 1) outputBuffer += ']';
			}


			//if(this->header.getSampleNumber() > 1) outputBuffer += ']';
		}
	}
}

void VariantReader::printFORMATCustomVector(buffer_type& outputBuffer,
					   const char& delimiter,
					   const U32& position,
					   const objects_type& objects,
					   std::vector<core::GTObject>& genotypes_unpermuted) const
{
	if(block_settings.format_all.display || this->variant_container.getBlock().n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta_container->at(position).format_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// Print key map
			outputBuffer += objects.format_field_names[targetKeys[0]];
			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ':';
				outputBuffer += objects.format_field_names[targetKeys[i]];
			};
			outputBuffer += delimiter;

			// First key
			// Cycle over keys
			objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
			// Cycle over samples
			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += ',';
				objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					outputBuffer += ',';
					objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
				}
			}
		}
	}
}

void VariantReader::printFORMATCustomVectorJSON(buffer_type& outputBuffer,
					   const char& delimiter,
					   const U32& position,
					   const objects_type& objects,
					   std::vector<core::GTObject>& genotypes_unpermuted) const
{
	if(block_settings.format_all.display || this->variant_container.getBlock().n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta_container->at(position).format_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// First key
			// Cycle over keys
			outputBuffer += "\"FORMAT-";
			outputBuffer += objects.format_field_names[targetKeys[0]];
			outputBuffer += '"';
			outputBuffer += ':';
			if(this->header.getSampleNumber() > 1) outputBuffer += '[';
			objects.format_containers[targetKeys[0]]->to_json_string(outputBuffer, position, 0);
			// Cycle over samples
			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += ',';
				objects.format_containers[targetKeys[0]]->to_json_string(outputBuffer, position, s);
			}
			if(this->header.getSampleNumber() > 1) outputBuffer += ']';

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"FORMAT-";
				outputBuffer += objects.format_field_names[targetKeys[i]];
				outputBuffer += '"';
				outputBuffer += ':';
				if(this->header.getSampleNumber() > 1) outputBuffer += '[';
				objects.format_containers[targetKeys[i]]->to_json_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					outputBuffer += ',';
					objects.format_containers[targetKeys[i]]->to_json_string(outputBuffer, position, s);
				}
				if(this->header.getSampleNumber() > 1) outputBuffer += ']';
			}
		}
	}
}

void VariantReader::printINFOVCF(buffer_type& outputBuffer,
				const char& delimiter,
				  const U32& position,
				  const objects_type& objects) const
{
	// Check if any INFO data exists at all
	if(this->variant_container.getBlock().footer.n_info_patterns == 0){
		return;
	}

	if(block_settings.info_all.display || this->variant_container.getBlock().n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta_container->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// First
			outputBuffer += objects.info_field_names[targetKeys[0]];
			if(objects.info_containers[targetKeys[0]]->emptyPosition(position) == false){
				outputBuffer += '=';
				objects.info_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ";";
				outputBuffer += objects.info_field_names[targetKeys[i]];
				if(this->header.info_fields[this->variant_container.getBlock().footer.info_offsets[targetKeys[i]].data_header.global_key].primitive_type == YON_VCF_HEADER_FLAG){
					continue;
				}
				if(objects.info_containers[targetKeys[i]]->emptyPosition(position)) continue;
				outputBuffer += '=';
				objects.info_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position);
			}
		} else {
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;
			outputBuffer += '.';
		}
	}
}

void VariantReader::printINFOVCF(buffer_type& outputBuffer,
				  const U32& position,
				  const objects_type& objects) const
{
	return(this->printINFOVCF(outputBuffer, '\t', position, objects));
}

void VariantReader::printINFOCustom(buffer_type& outputBuffer,
					 const char& delimiter,
					 const U32& position,
					 const objects_type& objects) const
{
	if(block_settings.info_all.display || this->variant_container.getBlock().n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta_container->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// Check if this target container is a FLAG
			if(this->header.info_fields[this->variant_container.getBlock().info_containers[targetKeys[0]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
				outputBuffer += objects.info_field_names[targetKeys[0]];
			} else {
				// Check if the positon is empty
				if(objects.info_containers[targetKeys[0]]->emptyPosition(position) == false){
					objects.info_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position);
				}
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				if(this->header.info_fields[this->variant_container.getBlock().info_containers[targetKeys[i]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
					outputBuffer +=objects.info_field_names[targetKeys[i]];
					continue;
				}
				if(objects.info_containers[targetKeys[i]]->emptyPosition(position)) continue;
				objects.info_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position);
			}
		}
	}
}

void VariantReader::printINFOCustomJSON(buffer_type& outputBuffer,
						 const char& delimiter,
						 const U32& position,
						 const objects_type& objects) const
{
	if(block_settings.info_all.display || this->variant_container.getBlock().n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta_container->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Check if this target container is a FLAG
			outputBuffer += "\"INFO-";
			outputBuffer += objects.info_field_names[targetKeys[0]];
			outputBuffer += "\":";
			if(this->header.info_fields[this->variant_container.getBlock().info_containers[targetKeys[0]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
				outputBuffer += "true";
			} else {
				objects.info_containers[targetKeys[0]]->to_json_string(outputBuffer, position);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"INFO-";
				outputBuffer += objects.info_field_names[targetKeys[i]];
				outputBuffer += "\":";
				if(this->header.info_fields[this->variant_container.getBlock().info_containers[targetKeys[i]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
					outputBuffer += "true";
					continue;
				}
				objects.info_containers[targetKeys[i]]->to_json_string(outputBuffer, position);
			}
		}
	}
}

TACHYON_VARIANT_CLASSIFICATION_TYPE VariantReader::classifyVariant(const meta_entry_type& meta, const U32& allele) const{
	const S32 ref_size = meta.alleles[0].size();
	const S32 diff = ref_size - meta.alleles[allele].size();
	//std::cerr << diff << ",";
	if(meta.alleles[0].allele[0] == '<' || meta.alleles[allele].allele[0] == '<') return(YON_VARIANT_CLASS_SV);
	else if(diff == 0){
		if(ref_size == 1 && meta.alleles[0].allele[0] != meta.alleles[allele].allele[0]){
			if(meta.alleles[allele].allele[0] == 'A' || meta.alleles[allele].allele[0] == 'T' || meta.alleles[allele].allele[0] == 'G' || meta.alleles[allele].allele[0] == 'C')
				return(YON_VARIANT_CLASS_SNP);
			else return(YON_VARIANT_CLASS_UNKNOWN);
		}
		else if(ref_size != 1){
			U32 characters_identical = 0;
			const U32 length_shortest = ref_size < meta.alleles[allele].size() ? ref_size : meta.alleles[allele].size();

			for(U32 c = 0; c < length_shortest; ++c){
				characters_identical += (meta.alleles[0].allele[c] == meta.alleles[allele].allele[c]);
			}

			if(characters_identical == 0) return(YON_VARIANT_CLASS_MNP);
			else return(YON_VARIANT_CLASS_CLUMPED);
		}
	} else {
		const U32 length_shortest = ref_size < meta.alleles[allele].size() ? ref_size : meta.alleles[allele].size();
		U32 characters_non_standard = 0;
		for(U32 c = 0; c < length_shortest; ++c){
			characters_non_standard += (meta.alleles[allele].allele[c] != 'A' && meta.alleles[allele].allele[c] != 'T' && meta.alleles[allele].allele[c] != 'C' && meta.alleles[allele].allele[c] !='G');
		}
		if(characters_non_standard) return(YON_VARIANT_CLASS_UNKNOWN);
		else return(YON_VARIANT_CLASS_INDEL);
	}
	return(YON_VARIANT_CLASS_UNKNOWN);
}

/**<
 * Outputs
 * @return
 */
const U64 VariantReader::outputVCF(void){
	U64 n_variants = 0;

	if(this->block_settings.annotate_extra){
		// fixme
		// if special
		// "FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"
		if(this->header.getInfoField("FS_A") == nullptr)          this->header.literals += "\n##INFO=<ID=FS_A,Number=A,Type=Float>";
		if(this->header.getInfoField("AN") == nullptr)            this->header.literals += "\n##INFO=<ID=AN,Number=A,Type=Integer>";
		if(this->header.getInfoField("NM") == nullptr)            this->header.literals += "\n##INFO=<ID=NM,Number=A,Type=Integer>";
		if(this->header.getInfoField("NPM") == nullptr)           this->header.literals += "\n##INFO=<ID=NPM,Number=A,Type=Integer>";
		if(this->header.getInfoField("AC") == nullptr)            this->header.literals += "\n##INFO=<ID=AC,Number=A,Type=Integer>";
		if(this->header.getInfoField("AC_FWD") == nullptr)        this->header.literals += "\n##INFO=<ID=AC_FWD,Number=A,Type=Integer>";
		if(this->header.getInfoField("AC_REV") == nullptr)        this->header.literals += "\n##INFO=<ID=AC_REV,Number=A,Type=Integer>";
		if(this->header.getInfoField("HWE_P") == nullptr)         this->header.literals += "\n##INFO=<ID=HWE_P,Number=A,Type=Float>";
		if(this->header.getInfoField("VT") == nullptr)            this->header.literals += "\n##INFO=<ID=VT,Number=A,Type=String>";
		if(this->header.getInfoField("AF") == nullptr)            this->header.literals += "\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">";
		if(this->header.getInfoField("MULTI_ALLELIC") == nullptr) this->header.literals += "\n##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag>";
		if(this->header.getInfoField("F_PIC") == nullptr)         this->header.literals += "\n##INFO=<ID=F_PIC,Number=A,Type=Float,Description=\"Population inbreeding coefficient (F-statistics)\">";
	}

	this->header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	this->header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

	this->header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + '\n';
	this->header.literals += this->getSettings().get_settings_string();

	// Output VCF header
	if(this->block_settings.show_vcf_header){
		this->header.writeVCFHeaderString(std::cout, this->block_settings.format_all.load || this->block_settings.format_list.size());
	}

	// If seek is active for targetted intervals
	if(this->interval_container.hasIntervals()){
		if(this->interval_container.build(this->header) == false)
			return false;

		if(this->interval_container.getBlockList().size()){
			for(U32 i = 0; i < this->interval_container.getBlockList().size(); ++i){
				if(this->getBlock(this->interval_container.getBlockList()[i]) == false){
					return(0);
				}
				n_variants += this->outputBlockVCF();
			}

			return(n_variants);
		} else { // Provided intervals but no matching YON blocks
			return(0);
		}
	}

	// While there are YON blocks
	while(this->nextBlock()) n_variants += this->outputBlockVCF();
	return(n_variants);
}

/**<
 *
 * @return
 */
const U64 VariantReader::outputCustom(void){
	U64 n_variants = 0;

	// While there are YON blocks
	while(this->nextBlock()) n_variants += this->outputBlockCustom();
	return(n_variants);
}

/**<
 *
 * @return
 */
const U32 VariantReader::outputBlockVCF(void){
	objects_type objects;
	this->loadObjects(objects);
	this->variant_filters.build();

	// Reserve memory for output buffer
	// This is much faster than writing directly to ostream because of synchronisation
	io::BasicBuffer output_buffer(256000);
	if(this->block_settings.format_all.load) output_buffer.resize(256000 + this->header.getSampleNumber()*2);

	// Todo: in cases of non-diploid
	std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());
	for(U32 i = 0; i < this->header.getSampleNumber(); ++i){
		genotypes_unpermuted[i].alleles = new core::GTObjectAllele;
	}

	// Print functionality
	print_format_function print_format = &self_type::printFORMATDummy;
	if(this->block_settings.format_ID_list.size()) print_format = &self_type::printFORMATCustom;
	else if(block_settings.format_all.display)     print_format = &self_type::printFORMATVCF;
	print_info_function   print_info   = &self_type::printINFOVCF;
	print_meta_function   print_meta   = &utility::to_vcf_string;
	print_filter_function print_filter = &self_type::printFILTER;

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::filterIntervalsDummy;
	if(this->interval_container.size()) filter_intervals = &self_type::filterIntervals;

	// Cycling over loaded meta objects
	for(U32 p = 0; p < objects.meta_container->size(); ++p){
		if(this->variant_filters.filter(objects, p) == false)
			continue;

		if((this->*filter_intervals)(objects.meta_container->at(p)) == false)
			continue;


		if(this->block_settings.custom_output_format)
			utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->header, this->block_settings);
		else
			utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->header);

		// Filter options
		(this->*print_filter)(output_buffer, p, objects);
		(this->*print_info)(output_buffer, '\t', p, objects);
		if(this->block_settings.annotate_extra) this->getGenotypeSummary(output_buffer, p, objects); // Todo: fixme
		(this->*print_format)(output_buffer, '\t', p, objects, genotypes_unpermuted);
		output_buffer += '\n';

		if(output_buffer.size() > 65536){
			std::cout.write(output_buffer.data(), output_buffer.size());
			output_buffer.reset();
			std::cout.flush();
		}
	}

	std::cout.write(output_buffer.data(), output_buffer.size());
	output_buffer.reset();
	std::cout.flush();

	return(objects.meta_container->size());
}

/**<
 *
 * @return
 */
const U32 VariantReader::outputBlockCustom(void) const{
	objects_type objects;
	this->loadObjects(objects);

	// Reserve memory for output buffer
	// This is much faster than writing directly to ostream because of syncing
	io::BasicBuffer output_buffer(256000 + this->header.getSampleNumber()*2);
	std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());

	// Todo: move to function
	//U32 info_match_limit = 1; // any match
	//info_match_limit = this->block_settings.info_list.size(); // all match

	// Function pointer to use
	print_format_function print_format = &self_type::printFORMATDummy;
	print_info_function   print_info   = &self_type::printINFODummy;
	print_meta_function   print_meta   = &utility::to_vcf_string;
	print_filter_function print_filter = &self_type::printFILTERDummy;
	filter_intervals_function filter_intervals = &self_type::filterIntervalsDummy;

	if(block_settings.output_json) print_meta = &utility::to_json_string;

	if(block_settings.format_all.display || this->variant_container.getBlock().n_format_loaded){
		if(block_settings.output_json){
			print_format = &self_type::printFORMATCustomVectorJSON;
		} else {
			if(block_settings.output_format_vector) print_format = &self_type::printFORMATCustomVector;
			else print_format = &self_type::printFORMATCustom;
		}
	}

	if(block_settings.info_all.display || this->variant_container.getBlock().n_info_loaded){
		if(block_settings.output_json) print_info = &self_type::printINFOCustomJSON;
		else print_info = &self_type::printINFOCustom;
	}

	if(block_settings.display_filter){
		if(block_settings.output_json) print_filter = &self_type::printFILTERJSON;
		else print_filter = &self_type::printFILTERCustom;
	}

	U32 n_records_returned = 0;

	if(block_settings.output_json) output_buffer += "\"block\":[";
	for(U32 position = 0; position < objects.meta_container->size(); ++position){
		//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
		//	continue;

		if(block_settings.output_json){
			if(position != 0) output_buffer += ",\n";
			output_buffer += "{";
		}
		++n_records_returned;

		(*print_meta)(output_buffer, this->block_settings.custom_delimiter_char, objects.meta_container->at(position), this->header, this->block_settings);
		(this->*print_filter)(output_buffer, position, objects);
		(this->*print_info)(output_buffer, this->block_settings.custom_delimiter_char, position, objects);
		(this->*print_format)(output_buffer, this->block_settings.custom_delimiter_char, position, objects, genotypes_unpermuted);

		if(block_settings.output_json) output_buffer += "}";
		else output_buffer += '\n';
		//output_buffer += "}";

		// Flush if buffer is large
		if(output_buffer.size() > 65536){
			std::cout.write(output_buffer.data(), output_buffer.size());
			output_buffer.reset();
			std::cout.flush();
		}
	}
	if(block_settings.output_json) output_buffer += "]";

	// Flush buffer
	std::cout.write(output_buffer.data(), output_buffer.size());
	output_buffer.reset();
	std::cout.flush();

	return(n_records_returned);
}

}
