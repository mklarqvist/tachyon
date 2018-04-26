#include "variant_reader.h"

namespace tachyon{

VariantReader::VariantReader() :
	filesize(0)
{}

VariantReader::VariantReader(const std::string& filename) :
	input_file(filename),
	filesize(0)
{}

VariantReader::~VariantReader(){}

VariantReader::VariantReader(const self_type& other) :
	input_file(other.input_file),
	filesize(other.filesize),
	settings(other.settings),
	header(other.header),
	footer(other.footer),
	index(other.index),
	checksums(other.checksums),
	keychain(other.keychain)
{
	this->stream.open(this->input_file, std::ios::in | std::ios::binary);
}

bool VariantReader::open(void){
	if(this->input_file.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	this->stream.open(this->input_file, std::ios::binary | std::ios::in | std::ios::ate);
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
	this->block.clear();

	if(!this->block.readHeaderFooter(this->stream))
		return false;

	if(!this->codec_manager.zstd_codec.decompress(this->block.footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	this->block.footer_support.buffer_data_uncompressed >> this->block.footer;
	this->parseSettings();

	// Attempts to read a YON block with the provided
	if(!this->block.read(this->stream, this->settings))
		return false;

	// encryption manager ascertainment
	if(this->block.header.controller.anyEncrypted){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption::EncryptionDecorator e;
		if(!e.decryptAES256(this->block, this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->block)){
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
	if(!block.read(this->stream, this->settings)){
		this->stream.seekg(position);
		return block;
	}

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->block)){
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
	objects.meta = new meta_container_type(this->block);
	if(this->block.header.controller.hasGT && settings.load_genotypes_all){
		objects.genotypes = new gt_container_type(this->block, *objects.meta);
	}

	objects.n_loaded_format = this->block.n_format_loaded;
	objects.format_fields = new format_interface_type*[objects.n_loaded_format];

	if(this->block.n_format_loaded){
		for(U32 i = 0; i < this->block.n_format_loaded; ++i){
			const U32 global_key = settings.load_format_ID_loaded[i].offset->data_header.global_key;
			std::vector<bool> matches = this->get_format_field_pattern_matches(this->header.format_fields[global_key].ID);

			if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.format_fields[i] = new containers::FormatContainer<S32>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header.format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.format_fields[i] = new containers::FormatContainer<std::string>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.format_fields[i] = new containers::FormatContainer<float>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			} else {
				objects.format_fields[i] = new containers::FormatContainer<U32>;
				objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
			}
		}
	}

	// Store as double pointers to avoid memory collisions because
	// info containers have different class members
	objects.n_loaded_info = this->block.n_info_loaded;
	objects.info_fields = new info_interface_type*[objects.n_loaded_info];

	if(this->block.n_info_loaded){
		for(U32 i = 0; i < this->block.n_info_loaded; ++i){
			const U32 global_key = settings.load_info_ID_loaded[i].offset->data_header.global_key;
			std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[global_key].ID);

			if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
				objects.info_fields[i] = new containers::InfoContainer<S32>(this->block.info_containers[i], *objects.meta, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
					  this->header.info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
				objects.info_fields[i] = new containers::InfoContainer<std::string>(this->block.info_containers[i], *objects.meta, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
				objects.info_fields[i] = new containers::InfoContainer<float>(this->block.info_containers[i], *objects.meta, matches);
				objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
			} else {
				objects.info_fields[i] = new containers::InfoContainer<U32>();
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
	objects.info_keep = std::vector<U32>(this->block.footer.n_info_patterns, 0);
	objects.format_keep = std::vector<U32>(this->block.footer.n_format_patterns, 0);

	// This vector of vectors keeps the local INFO identifiers for the matched global
	// identifiers in a given hash pattern
	//
	// For example: x[5] contains the local IDs for loaded INFO streams for pattern ID 5
	objects.local_match_keychain_info = std::vector< std::vector<U32> >(this->block.footer.n_info_patterns);
	objects.local_match_keychain_format = std::vector< std::vector<U32> >(this->block.footer.n_format_patterns);

	// If loading all INFO values then return them in the ORIGINAL order
	if(this->settings.load_info){
		for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
					if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_keep[i];
					}
				}
			}
		}
	}
	// If loading custom INFO fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_info[i].push_back(k);
						++objects.info_keep[i];
					}
				}
			}
		}
	}

	// For FORMAT
	// If loading all FORMAT values then return them in the ORIGINAL order
	if(this->settings.load_format){
		for(U32 i = 0; i < this->block.footer.n_format_patterns; ++i){ // Number of info patterns
			for(U32 j = 0; j < this->block.footer.format_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
				for(U32 k = 0; k < settings.load_format_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
					if(this->block.footer.format_offsets[this->block.footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_format_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_format[i].push_back(k);
						++objects.format_keep[i];
					}
				}
			}
		}
	}
	// If loading custom FORMAT fields then return them in the REQUESTED order
	else {
		for(U32 i = 0; i < this->block.footer.n_format_patterns; ++i){ // i = Number of info patterns
			for(U32 k = 0; k < settings.load_format_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
				for(U32 j = 0; j < this->block.footer.format_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
					if(this->block.footer.format_offsets[this->block.footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_format_ID_loaded[k].offset->data_header.global_key){
						objects.local_match_keychain_format[i].push_back(k);
						++objects.format_keep[i];
					}
				}
			}
		}
	}

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
		for(U32 i = 0; i < this->block.footer.n_info_streams; ++i){
			if(this->block.footer.info_offsets[i].data_header.global_key == global_info_value){
				local_key = i;
			}
		}

		if(local_key == -1){
			std::cerr << "could not find local" << std::endl;
			return ret;
		}

		ret.resize(this->block.footer.n_info_patterns, false);
		for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){
			//std::cerr << i << '\t' << this->block.footer.info_bit_vectors[i][local_key] << std::endl;
			ret[i] = this->block.footer.info_bit_vectors[i][local_key];
		}
	}
	return(ret);
}

const std::vector<bool> VariantReader::get_format_field_pattern_matches(const std::string& field_name) const{
	int global_format_value = this->has_format_field(field_name);
	std::vector<bool> ret;
	if(global_format_value >= 0){
		S32 local_key = -1;
		for(U32 i = 0; i < this->block.footer.n_format_streams; ++i){
			if(this->block.footer.format_offsets[i].data_header.global_key == global_format_value){
				local_key = i;
			}
		}

		if(local_key == -1){
			std::cerr << "could not find local" << std::endl;
			return ret;
		}

		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		ret.resize(this->block.footer.n_format_patterns, false);
		for(U32 i = 0; i < this->block.footer.n_format_patterns; ++i){
			//std::cerr << i << '\t' << this->block.index_entry.format_bit_vectors[i][local_format_field_id] << std::endl;
			ret[i] = this->block.footer.format_bit_vectors[i][local_key];
		}
	}
	return(ret);
}


void VariantReader::printFILTER(buffer_type& outputBuffer,
				 const U32& position,
				 const objects_type& objects) const
{
	if(this->block.footer.n_filter_streams){
		const U32& n_filter_keys = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			if(outputBuffer.back() != '\t') outputBuffer += '\t';

			// Local key -> global key
			outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
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
	if(this->block.footer.n_filter_streams){
		const U32& n_filter_keys = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			if(outputBuffer.back() != this->settings.custom_delimiter_char) outputBuffer += this->settings.custom_delimiter_char;

			// Local key -> global key
			outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
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
	if(this->block.footer.n_filter_streams){
		const U32& n_filter_keys = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].n_keys;
		const U32* filter_keys   = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].local_keys;
		if(n_filter_keys){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Local key -> global key
			outputBuffer += "\"FILTER-";
			outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
			outputBuffer += "\":true";
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ',';
				outputBuffer += "\"FILTER-";
				outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
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
	if(settings.load_format && this->block.n_format_loaded){
		if(this->block.n_format_loaded){
			const U32& n_format_keys = this->block.footer.format_bit_vectors[objects.meta->at(position).format_pattern_id].n_keys;
			const U32* format_keys   = this->block.footer.format_bit_vectors[objects.meta->at(position).format_pattern_id].local_keys;
			if(n_format_keys){
				// Print key map
				buffer += this->header.format_fields[this->block.footer.format_offsets[format_keys[0]].data_header.global_key].ID;
				for(U32 i = 1; i < n_format_keys; ++i){
					buffer += ':';
					buffer += this->header.format_fields[this->block.footer.format_offsets[format_keys[i]].data_header.global_key].ID;
				}
				buffer += delimiter;

				// Todo: print if no GT data
				// Begin print FORMAT data for each sample
				if(this->settings.load_ppa && this->block.header.controller.hasGTPermuted){
					objects.genotypes->at(position).getObjects(genotypes_unpermuted, this->header.getSampleNumber(), this->block.ppa_manager);
				} else {
					objects.genotypes->at(position).getObjects(genotypes_unpermuted, this->header.getSampleNumber());
				}

				buffer << genotypes_unpermuted[0];
				for(U32 i = 1; i < n_format_keys; ++i){
					buffer += ':';
					objects.format_fields[format_keys[i]]->to_vcf_string(buffer, position, 0);
				}

				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					buffer += delimiter;
					buffer << genotypes_unpermuted[s];
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer  += ':';
						objects.format_fields[format_keys[i]]->to_vcf_string(buffer, position, s);
					}
				}
			} else { // have no keys
				buffer += ".";
				buffer += delimiter;
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
	if(settings.load_format || this->block.n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta->at(position).format_pattern_id];
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
			objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ':';
				objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
			}
			//if(targetKeys.size() > 1) outputBuffer += ']';

			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += delimiter;
				//if(targetKeys.size() > 1) outputBuffer += '[';
				objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer  += ':';
					objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
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
	if(settings.load_format || this->block.n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta->at(position).format_pattern_id];
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
			objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
			// Cycle over samples
			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += ',';
				objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					outputBuffer += ',';
					objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
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
	if(settings.load_format || this->block.n_format_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta->at(position).format_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// First key
			// Cycle over keys
			outputBuffer += "\"FORMAT-";
			outputBuffer += objects.format_field_names[targetKeys[0]];
			outputBuffer += '"';
			outputBuffer += ':';
			if(this->header.getSampleNumber() > 1) outputBuffer += '[';
			objects.format_fields[targetKeys[0]]->to_json_string(outputBuffer, position, 0);
			// Cycle over samples
			for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
				outputBuffer += ',';
				objects.format_fields[targetKeys[0]]->to_json_string(outputBuffer, position, s);
			}
			if(this->header.getSampleNumber() > 1) outputBuffer += ']';

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"FORMAT-";
				outputBuffer += objects.format_field_names[targetKeys[i]];
				outputBuffer += '"';
				outputBuffer += ':';
				if(this->header.getSampleNumber() > 1) outputBuffer += '[';
				objects.format_fields[targetKeys[i]]->to_json_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					outputBuffer += ',';
					objects.format_fields[targetKeys[i]]->to_json_string(outputBuffer, position, s);
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
	if(this->block.footer.n_info_patterns == 0){
		outputBuffer += '.';
		return;
	}

	if(settings.load_info || this->block.n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta->at(position).info_pattern_id];
		if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

		if(targetKeys.size()){
			// First
			outputBuffer += objects.info_field_names[targetKeys[0]];
			if(objects.info_fields[targetKeys[0]]->emptyPosition(position) == false){
				outputBuffer += '=';
				objects.info_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ";";
				outputBuffer += objects.info_field_names[targetKeys[i]];
				if(this->header.info_fields[this->block.footer.info_offsets[targetKeys[i]].data_header.global_key].primitive_type == YON_VCF_HEADER_FLAG){
					continue;
				}
				if(objects.info_fields[targetKeys[i]]->emptyPosition(position)) continue;
				outputBuffer += '=';
				objects.info_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position);
			}
		} else {
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
	if(settings.load_info || this->block.n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// Check if this target container is a FLAG
			if(this->header.info_fields[this->block.info_containers[targetKeys[0]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
				outputBuffer += objects.info_field_names[targetKeys[0]];
			} else {
				// Check if the positon is empty
				if(objects.info_fields[targetKeys[0]]->emptyPosition(position) == false){
					objects.info_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position);
				}
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				if(this->header.info_fields[this->block.info_containers[targetKeys[i]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
					outputBuffer +=objects.info_field_names[targetKeys[i]];
					continue;
				}
				if(objects.info_fields[targetKeys[i]]->emptyPosition(position)) continue;
				objects.info_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position);
			}
		}
	}
}

void VariantReader::printINFOCustomJSON(buffer_type& outputBuffer,
						 const char& delimiter,
						 const U32& position,
						 const objects_type& objects) const
{
	if(settings.load_info || this->block.n_info_loaded){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Check if this target container is a FLAG
			outputBuffer += "\"INFO-";
			outputBuffer += objects.info_field_names[targetKeys[0]];
			outputBuffer += "\":";
			if(this->header.info_fields[this->block.info_containers[targetKeys[0]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
				outputBuffer += "true";
			} else {
				objects.info_fields[targetKeys[0]]->to_json_string(outputBuffer, position);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"INFO-";
				outputBuffer += objects.info_field_names[targetKeys[i]];
				outputBuffer += "\":";
				if(this->header.info_fields[this->block.info_containers[targetKeys[i]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
					outputBuffer += "true";
					continue;
				}
				objects.info_fields[targetKeys[i]]->to_json_string(outputBuffer, position);
			}
		}
	}
}

}
