#include "variant_reader.h"

namespace tachyon{

VariantReader::VariantReader() :
	filesize(0)
{}

VariantReader::VariantReader(const std::string& filename) :
	input_file(filename),
	filesize(0)
{}

// Dtor
VariantReader::~VariantReader(){}

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
	this->stream >> this->header;
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
		std::cerr << utility::timestamp("ERROR") << "Corrupted!" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->stream.tellg() == this->footer.offset_end_of_data)
		return false;

	// Reset and re-use
	this->block.clear();

	// Attempts to read a YON block with the provided
	// settings
	if(!this->block.read(this->stream, this->settings))
		return false;

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->block))
		return false;

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
			//std::cerr << i << '\t' << this->block.index_entry.info_bit_vectors[i][local_info_field_id] << std::endl;
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

}
