#include "variant_reader.h"

namespace tachyon{

VariantReader::VariantReader()
{}

VariantReader::VariantReader(const std::string& filename) :
	basic_reader(filename)
{}

VariantReader::~VariantReader(){}

VariantReader::VariantReader(const self_type& other) :
	basic_reader(other.basic_reader),
	block_settings(other.block_settings),
	settings(other.settings),
	global_header(other.global_header),
	global_footer(other.global_footer),
	index(other.index),
	checksums(other.checksums),
	keychain(other.keychain)
{
	this->basic_reader.open();
}

bool VariantReader::open(void){
	if(this->settings.input.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	if(this->basic_reader.open() == false){
		std::cerr << "Failed to open" << std::endl;
		return false;
	}

	//this->basic_reader.stream_.open(this->settings.input, std::ios::binary | std::ios::in | std::ios::ate);
	//this->filesize = (U64)this->basic_reader.stream_.tellg();

	if(this->basic_reader.filesize_ <= YON_FOOTER_LENGTH){
		std::cerr << utility::timestamp("ERROR") << "File is corrupted!" << std::endl;
		return false;
	}

	// Seek to start of footer
	this->basic_reader.stream_.seekg((U64)this->basic_reader.filesize_ - YON_FOOTER_LENGTH);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to seek in file!" << std::endl;
		return false;
	}
	this->basic_reader.stream_ >> this->global_footer;

	// Validate footer
	if(this->global_footer.validate() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate footer!" << std::endl;
		return false;
	}

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to read file!" << std::endl;
		return false;
	}

	// Seek to start of file
	this->basic_reader.stream_.seekg(0);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to rewind file!" << std::endl;
		return false;
	}

	// Load header
	//this->stream >> this->global_header;
	char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	this->basic_reader.stream_.read(&magic_string[0], tachyon::constants::FILE_HEADER_LENGTH);
	if(strncmp(&magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate Tachyon magic string!" << std::endl;
		return false;
	}

	uint32_t l_data   = 0;
	uint32_t l_c_data = 0;
	utility::DeserializePrimitive(l_data, this->basic_reader.stream_);
	utility::DeserializePrimitive(l_c_data, this->basic_reader.stream_);

	io::BasicBuffer header_uncompressed(l_data + 1024);
	io::BasicBuffer header_compressed(l_c_data + 1024); header_compressed.n_chars   = l_c_data;

	this->basic_reader.stream_.read(header_compressed.data(), l_c_data);

	if(!this->codec_manager.zstd_codec.Decompress(header_compressed, header_uncompressed)){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	assert(header_uncompressed.size() == l_data);
	header_uncompressed >> this->global_header; // parse header from buffer

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to get header!" << std::endl;
		return false;
	}

	this->variant_container << this->global_header;

	// Keep track of start position
	const U64 return_pos = this->basic_reader.stream_.tellg();
	this->basic_reader.stream_.seekg(this->global_footer.offset_end_of_data);
	this->basic_reader.stream_ >> this->index;
	this->basic_reader.stream_ >> this->checksums;
	this->basic_reader.stream_.seekg(return_pos);

	// Parse settings
	this->getBlockSettings().parseSettings(this->global_header);

	return(this->basic_reader.stream_.good());
}

bool VariantReader::nextBlock(){
	// If the stream is faulty then return
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->basic_reader.stream_.tellg() == this->global_footer.offset_end_of_data)
		return false;

	// Reset and re-use
	this->variant_container.reset();

	if(!this->variant_container.getBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return false;


	if(!this->codec_manager.zstd_codec.Decompress(this->variant_container.getBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	this->variant_container.getBlock().footer_support.buffer_data_uncompressed >> this->variant_container.getBlock().footer;

	// Attempts to read a YON block with the settings provided
	if(!this->variant_container.readBlock(this->basic_reader.stream_, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.anyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.decryptAES256(this->variant_container.getBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}
	//assert(this->variant_container.getBlock().gt_ppa != nullptr);
	this->variant_container.getBlock().gt_ppa = new yon_gt_ppa;
	this->variant_container.getBlock().gt_ppa->n_samples = this->global_header.GetNumberSamples();

	// Internally decompress available data
	if(!this->codec_manager.Decompress(this->variant_container.getBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

bool VariantReader::getBlock(const index_entry_type& index_entry){
	// If the stream is faulty then return
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// Reset and re-use
	this->variant_container.getBlock().clear();

	// Seek to target block ID with the help of the index
	this->basic_reader.stream_.seekg(index_entry.byte_offset);
	if(this->basic_reader.stream_.good() == false){
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to seek to given offset using target index entry!" << std::endl;
		return(false);
	}

	if(!this->variant_container.getBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return false;

	if(!this->codec_manager.zstd_codec.Decompress(this->variant_container.getBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
		return(false);
	}
	this->variant_container.getBlock().footer_support.buffer_data_uncompressed >> this->variant_container.getBlock().footer;

	// Attempts to read a YON block with the provided
	if(!this->variant_container.readBlock(this->basic_reader.stream_, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.getBlock().header.controller.anyEncrypted){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.decryptAES256(this->variant_container.getBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.Decompress(this->variant_container.getBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

containers::VariantBlockContainer VariantReader::getBlock(){
	// If the stream is faulty then return
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return variant_container_type();
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->basic_reader.stream_.tellg() == this->global_footer.offset_end_of_data)
		return variant_container_type();

	// Reset and re-use
	containers::VariantBlockContainer block;

	if(!block.getBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return variant_container_type();

	if(!this->codec_manager.zstd_codec.Decompress(block.getBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	block.getBlock().footer_support.buffer_data_uncompressed >> block.getBlock().footer;

	// Attempts to read a YON block with the settings provided
	if(!block.readBlock(this->basic_reader.stream_, this->getBlockSettings()))
		return variant_container_type();

	// encryption manager ascertainment
	if(block.anyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return variant_container_type();
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.decryptAES256(block.getBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return variant_container_type();
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.Decompress(block.getBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return variant_container_type();
	}

	// All passed
	return block;
}

void VariantReader::printFILTER(buffer_type& outputBuffer,
				 const U32& position,
				 const objects_type& objects) const
{
	if(outputBuffer.back() != '\t') outputBuffer += '\t';

	if(this->block_settings.display_filter && this->variant_container.getBlock().footer.n_filter_streams){
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.size();
		const int* filter_keys   = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.data();
		if(n_filter_keys){
			// Local key -> global key
			outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].id;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].id;
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
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.size();
		const int* filter_keys   = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.data();
		if(n_filter_keys){
			if(outputBuffer.back() != this->block_settings.custom_delimiter_char) outputBuffer += this->block_settings.custom_delimiter_char;

			// Local key -> global key
			outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].id;
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ';';
				outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].id;
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
		const U32& n_filter_keys = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.size();
		const int* filter_keys   = this->variant_container.getBlock().footer.filter_patterns[(*objects.meta_container)[position].filter_pattern_id].pattern.data();
		if(n_filter_keys){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Local key -> global key
			outputBuffer += "\"FILTER-";
			outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[0]].data_header.global_key].id;
			outputBuffer += "\":true";
			for(U32 i = 1; i < n_filter_keys; ++i){
				outputBuffer += ',';
				outputBuffer += "\"FILTER-";
				outputBuffer += this->global_header.filter_fields_[this->variant_container.getBlock().footer.filter_offsets[filter_keys[i]].data_header.global_key].id;
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
	if((this->block_settings.display_static & YON_BLK_BV_FORMAT) && objects.n_loaded_format){
		if(objects.n_loaded_format){
			const U32& n_format_keys = this->variant_container.getBlock().footer.format_patterns[objects.meta_container->at(position).format_pattern_id].pattern.size();
			const int* format_keys   = this->variant_container.getBlock().footer.format_patterns[objects.meta_container->at(position).format_pattern_id].pattern.data();
			if(n_format_keys){
				if(buffer.back() != delimiter) buffer += delimiter;

				// Print key map
				buffer += this->global_header.format_fields_[this->variant_container.getBlock().footer.format_offsets[format_keys[0]].data_header.global_key].id;
				for(U32 i = 1; i < n_format_keys; ++i){
					buffer += ':';
					buffer += this->global_header.format_fields_[this->variant_container.getBlock().footer.format_offsets[format_keys[i]].data_header.global_key].id;
				}
				buffer += delimiter;

				// Todo: print if no GT data
				// Begin print FORMAT data for each sample
				if(this->variant_container.getBlock().header.controller.hasGT){
					if((this->block_settings.display_static & YON_BLK_BV_PPA) && this->variant_container.getBlock().header.controller.hasGTPermuted){
						//objects.genotype_container->at(position).getObjects(this->global_header.GetNumberSamples(), genotypes_unpermuted, this->variant_container.getBlock().ppa_manager);
					} else {
						//objects.genotype_container->at(position).getObjects(this->global_header.GetNumberSamples(), genotypes_unpermuted);
					}

					buffer << genotypes_unpermuted[0];
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						objects.format_containers[format_keys[i]]->to_vcf_string(buffer, position, 0);
					}

					for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
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

					for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
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
	if((this->block_settings.display_static & YON_BLK_BV_FORMAT) || objects.n_loaded_format){
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
			//if(this->global_header.GetNumberSamples() > 1) outputBuffer += '[';

			//if(targetKeys.size() > 1) outputBuffer += '[';
			objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ':';
				objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
			}
			//if(targetKeys.size() > 1) outputBuffer += ']';

			for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
				outputBuffer += delimiter;
				//if(targetKeys.size() > 1) outputBuffer += '[';
				objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer  += ':';
					objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
				}
				//if(targetKeys.size() > 1) outputBuffer += ']';
			}


			//if(this->global_header.GetNumberSamples() > 1) outputBuffer += ']';
		}
	}
}

void VariantReader::printFORMATCustomVector(buffer_type& outputBuffer,
					   const char& delimiter,
					   const U32& position,
					   const objects_type& objects,
					   std::vector<core::GTObject>& genotypes_unpermuted) const
{
	if((this->block_settings.display_static & YON_BLK_BV_FORMAT) || objects.n_loaded_format){
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
			for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
				outputBuffer += ',';
				objects.format_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				objects.format_containers[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
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
	if((this->block_settings.display_static & YON_BLK_BV_FORMAT) || objects.n_loaded_format){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta_container->at(position).format_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// First key
			// Cycle over keys
			outputBuffer += "\"FORMAT-";
			outputBuffer += objects.format_field_names[targetKeys[0]];
			outputBuffer += '"';
			outputBuffer += ':';
			if(this->global_header.GetNumberSamples() > 1) outputBuffer += '[';
			objects.format_containers[targetKeys[0]]->to_json_string(outputBuffer, position, 0);
			// Cycle over samples
			for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
				outputBuffer += ',';
				objects.format_containers[targetKeys[0]]->to_json_string(outputBuffer, position, s);
			}
			if(this->global_header.GetNumberSamples() > 1) outputBuffer += ']';

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"FORMAT-";
				outputBuffer += objects.format_field_names[targetKeys[i]];
				outputBuffer += '"';
				outputBuffer += ':';
				if(this->global_header.GetNumberSamples() > 1) outputBuffer += '[';
				objects.format_containers[targetKeys[i]]->to_json_string(outputBuffer, position, 0);
				// Cycle over samples
				for(U64 s = 1; s < this->global_header.GetNumberSamples(); ++s){
					outputBuffer += ',';
					objects.format_containers[targetKeys[i]]->to_json_string(outputBuffer, position, s);
				}
				if(this->global_header.GetNumberSamples() > 1) outputBuffer += ']';
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
	if(this->variant_container.getBlock().footer.n_info_patterns == 0)
		return;

	if((this->block_settings.display_static & YON_BLK_BV_INFO) || objects.n_loaded_info){
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
				if(this->global_header.info_fields_[this->variant_container.getBlock().footer.info_offsets[targetKeys[i]].data_header.global_key].yon_type == YON_VCF_HEADER_FLAG)
					continue;

				if(objects.info_containers[targetKeys[i]]->emptyPosition(position))
					continue;

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
	if((this->block_settings.display_static & YON_BLK_BV_INFO) || objects.n_loaded_info){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta_container->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != delimiter) outputBuffer += delimiter;

			// Check if this target container is a FLAG
			if(this->global_header.info_fields_[this->variant_container.getBlock().info_containers[targetKeys[0]].header.getGlobalKey()].yon_type == YON_VCF_HEADER_FLAG){
				outputBuffer += objects.info_field_names[targetKeys[0]];
			} else {
				// Check if the position is empty
				if(objects.info_containers[targetKeys[0]]->emptyPosition(position) == false){
					objects.info_containers[targetKeys[0]]->to_vcf_string(outputBuffer, position);
				}
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += delimiter;
				if(this->global_header.info_fields_[this->variant_container.getBlock().info_containers[targetKeys[i]].header.getGlobalKey()].yon_type == YON_VCF_HEADER_FLAG){
					outputBuffer += objects.info_field_names[targetKeys[i]];
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
	if((this->block_settings.display_static & YON_BLK_BV_INFO) || objects.n_loaded_info){
		const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta_container->at(position).info_pattern_id];
		if(targetKeys.size()){
			if(outputBuffer.back() != ',') outputBuffer += ',';
			// Check if this target container is a FLAG
			outputBuffer += "\"INFO-";
			outputBuffer += objects.info_field_names[targetKeys[0]];
			outputBuffer += "\":";
			if(this->global_header.info_fields_[this->variant_container.getBlock().info_containers[targetKeys[0]].header.getGlobalKey()].yon_type == YON_VCF_HEADER_FLAG){
				outputBuffer += "true";
			} else {
				objects.info_containers[targetKeys[0]]->to_json_string(outputBuffer, position);
			}

			for(U32 i = 1; i < targetKeys.size(); ++i){
				outputBuffer += ',';
				outputBuffer += "\"INFO-";
				outputBuffer += objects.info_field_names[targetKeys[i]];
				outputBuffer += "\":";
				if(this->global_header.info_fields_[this->variant_container.getBlock().info_containers[targetKeys[i]].header.getGlobalKey()].yon_type == YON_VCF_HEADER_FLAG){
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

U64 VariantReader::OutputVcf(void){
	while(this->nextBlock()){
		objects_type& objects = *this->getCurrentBlock().loadObjects(this->block_settings);
		containers::yon1_t* entries = this->getCurrentBlock().LazyEvaluate(objects);

		// Todo: in cases of non-diploid
		io::BasicBuffer output_buffer(100000);
		/*
		std::vector<core::GTObject> genotypes_unpermuted(this->global_header.GetNumberSamples());
		for(U32 i = 0; i < this->global_header.GetNumberSamples(); ++i){
			genotypes_unpermuted[i].alleles = new core::GTObjectAllele;
		}

		for(U32 i = 0; i < objects.meta_container->size(); ++i){
			entries[i].gt->getObjects(this->global_header.GetNumberSamples(), genotypes_unpermuted);
			for(U32 j = 0; j < genotypes_unpermuted.size(); ++j){
				output_buffer << genotypes_unpermuted[j];
				output_buffer += ' ';
			}
			output_buffer += '\n';

			std::cout.write(output_buffer.data(), output_buffer.size());
			output_buffer.reset();
		}

		// temp
		delete [] entries;
		continue;
		*/

		for(U32 i = 0; i < objects.meta_container->size(); ++i){
			utility::to_vcf_string(output_buffer, '\t', *entries[i].meta, this->global_header);
			output_buffer += '\t';

			// Print Filter.
			const uint32_t n_filter_avail = entries[i].filter_ids->size();
			if(n_filter_avail){
				output_buffer += entries[i].filter_hdr[0]->id;
				for(U32 j = 1; j < n_filter_avail; ++j){
					output_buffer += ';';
					output_buffer += entries[i].filter_hdr[j]->id;
				}
			} else {
				output_buffer += '.';
			}
			output_buffer += '\t';

			// Print Info.
			const uint32_t n_info_avail = entries[i].info_ids->size();
			if(n_info_avail){
				if(entries[i].info_hdr[0]->yon_type == YON_VCF_HEADER_FLAG){
					output_buffer += entries[i].info_hdr[0]->id;
				} else {
					output_buffer += entries[i].info_hdr[0]->id;
					output_buffer += '=';
					entries[i].info[0]->to_vcf_string(output_buffer);
				}

				for(U32 j = 1; j < n_info_avail; ++j){
					output_buffer += ';';
					if(entries[i].info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
						output_buffer += entries[i].info_hdr[j]->id;
					} else {
						output_buffer += entries[i].info_hdr[j]->id;
						output_buffer += '=';
						entries[i].info[j]->to_vcf_string(output_buffer);
					}
				}
			}
			output_buffer += '\t';

			// Print Format.
			const uint32_t n_format_avail = entries[i].format_ids->size();
			if(n_format_avail){
				output_buffer += entries[i].format_hdr[0]->id;
				for(U32 j = 1; j < n_format_avail; ++j){
					output_buffer += ':';
					output_buffer += entries[i].format_hdr[j]->id;
				}
				output_buffer += '\t';

				// Todo: this is currently only valid if GT is available for
				// this record.

				if(n_format_avail > entries[i].is_loaded_gt){
					for(U32 s = 0; s < this->global_header.GetNumberSamples(); ++s){
						entries[i].format_containers[entries[i].is_loaded_gt]->to_vcf_string(output_buffer, i, s);
						for(U32 g = entries[i].is_loaded_gt + 1; g < n_format_avail; ++g){
							output_buffer += ':';
							entries[i].format_containers[g]->to_vcf_string(output_buffer, i, s);
						}
						output_buffer += '\t';
					}
				}
			}

			output_buffer += '\n';

			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
			}
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		delete [] entries;
	}

	return 0;
}

/**<
 * Outputs
 * @return
 */
U64 VariantReader::outputVCF(void){
	U64 n_variants = 0;

	if(this->block_settings.annotate_extra){
		// fixme
		// if special
		// "FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"
		if(this->global_header.GetInfo("FS_A") == nullptr)          this->global_header.literals_ += "\n##INFO=<ID=FS_A,Number=A,Type=Float>";
		if(this->global_header.GetInfo("AN") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AN,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("NM") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=NM,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("NPM") == nullptr)           this->global_header.literals_ += "\n##INFO=<ID=NPM,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AC,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC_FWD") == nullptr)        this->global_header.literals_ += "\n##INFO=<ID=AC_FWD,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC_REV") == nullptr)        this->global_header.literals_ += "\n##INFO=<ID=AC_REV,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("HWE_P") == nullptr)         this->global_header.literals_ += "\n##INFO=<ID=HWE_P,Number=A,Type=Float>";
		if(this->global_header.GetInfo("VT") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=VT,Number=A,Type=String>";
		if(this->global_header.GetInfo("AF") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">";
		if(this->global_header.GetInfo("MULTI_ALLELIC") == nullptr) this->global_header.literals_ += "\n##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag>";
		if(this->global_header.GetInfo("F_PIC") == nullptr)         this->global_header.literals_ += "\n##INFO=<ID=F_PIC,Number=A,Type=Float,Description=\"Population inbreeding coefficient (F-statistics)\">";
	}

	this->global_header.literals_ += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	this->global_header.literals_ += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

	this->global_header.literals_ += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + '\n';
	this->global_header.literals_ += this->getSettings().get_settings_string();

	// Output VCF header
	if(this->block_settings.show_vcf_header){
		//this->global_header.writeHeaderVCF(std::cout, this->block_settings.format_all.load || this->block_settings.format_list.size());
	}

	// If seek is active for targetted intervals
	if(this->interval_container.hasIntervals()){
		if(this->interval_container.build(this->global_header) == false)
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
U64 VariantReader::outputCustom(void){
	U64 n_variants = 0;

	// While there are YON blocks
	while(this->nextBlock()) n_variants += this->outputBlockCustom();
	return(n_variants);
}

/**<
 *
 * @return
 */
U32 VariantReader::outputBlockVCF(void){
	objects_type& objects = *this->getCurrentBlock().loadObjects(this->block_settings);

	// Reserve memory for output buffer
	// This is much faster than writing directly to ostream because of synchronisation
	io::BasicBuffer output_buffer(256000);
	if((this->block_settings.display_static & YON_BLK_BV_FORMAT)) output_buffer.resize(256000 + this->global_header.GetNumberSamples()*2);

	// Todo: in cases of non-diploid
	std::vector<core::GTObject> genotypes_unpermuted(this->global_header.GetNumberSamples());
	for(U32 i = 0; i < this->global_header.GetNumberSamples(); ++i){
		genotypes_unpermuted[i].alleles = new core::GTObjectAllele;
	}

	// Print functionality
	print_format_function print_format = &self_type::printFORMATDummy;
	if(this->block_settings.format_ID_list.size()) print_format = &self_type::printFORMATCustom;
	else if((this->block_settings.display_static & YON_BLK_BV_FORMAT))     print_format = &self_type::printFORMATVCF;
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
			utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->global_header, this->block_settings);
		else
			utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->global_header);

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
U32 VariantReader::outputBlockCustom(void){
	objects_type& objects = *this->getCurrentBlock().loadObjects(this->block_settings);

	// Reserve memory for output buffer
	// This is much faster than writing directly to ostream because of syncing
	io::BasicBuffer output_buffer(256000 + this->global_header.GetNumberSamples()*2);
	std::vector<core::GTObject> genotypes_unpermuted(this->global_header.GetNumberSamples());

	// Todo: move to function
	//U32 info_match_limit = 1; // any match
	//info_match_limit = this->block_settings.info_list.size(); // all match

	// Function pointer to use
	print_format_function print_format = &self_type::printFORMATDummy;
	print_info_function   print_info   = &self_type::printINFODummy;
	print_meta_function   print_meta   = &utility::to_vcf_string;
	print_filter_function print_filter = &self_type::printFILTERDummy;
	if(this->block_settings.output_json) print_meta = &utility::to_json_string;

	if((this->block_settings.display_static & YON_BLK_BV_FORMAT) || objects.n_loaded_format){
		if(this->block_settings.output_json){
			print_format = &self_type::printFORMATCustomVectorJSON;
		} else {
			if(this->block_settings.output_format_vector) print_format = &self_type::printFORMATCustomVector;
			else print_format = &self_type::printFORMATCustom;
		}
	}

	if((this->block_settings.display_static & YON_BLK_BV_INFO) || objects.n_loaded_info){
		if(this->block_settings.output_json) print_info = &self_type::printINFOCustomJSON;
		else print_info = &self_type::printINFOCustom;
	}

	if(this->block_settings.display_filter){
		if(this->block_settings.output_json) print_filter = &self_type::printFILTERJSON;
		else print_filter = &self_type::printFILTERCustom;
	}

	U32 n_records_returned = 0;

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::filterIntervalsDummy;
	if(this->interval_container.size()) filter_intervals = &self_type::filterIntervals;


	if(this->block_settings.output_json) output_buffer += "\"block\":[";
	for(U32 position = 0; position < objects.meta_container->size(); ++position){
		if(this->variant_filters.filter(objects, position) == false)
			continue;

		if((this->*filter_intervals)(objects.meta_container->at(position)) == false)
			continue;

		//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
		//	continue;

		if(this->block_settings.output_json){
			if(position != 0) output_buffer += ",\n";
			output_buffer += "{";
		}
		++n_records_returned;

		(*print_meta)(output_buffer, this->block_settings.custom_delimiter_char, objects.meta_container->at(position), this->global_header, this->block_settings);
		(this->*print_filter)(output_buffer, position, objects);
		(this->*print_info)(output_buffer, this->block_settings.custom_delimiter_char, position, objects);
		(this->*print_format)(output_buffer, this->block_settings.custom_delimiter_char, position, objects, genotypes_unpermuted);

		if(this->block_settings.output_json) output_buffer += "}";
		else output_buffer += '\n';
		//output_buffer += "}";

		// Flush if buffer is large
		if(output_buffer.size() > 65536){
			std::cout.write(output_buffer.data(), output_buffer.size());
			output_buffer.reset();
			std::cout.flush();
		}
	}
	if(this->block_settings.output_json) output_buffer += "]";

	// Flush buffer
	std::cout.write(output_buffer.data(), output_buffer.size());
	output_buffer.reset();
	std::cout.flush();

	return(n_records_returned);
}

}
