#include "tachyon_reader.h"

namespace tachyon{

TachyonReader::TachyonReader() :
	filesize(0)
{}

TachyonReader::TachyonReader(const std::string& filename) :
	input_file(filename),
	filesize(0)
{}

// Dtor
TachyonReader::~TachyonReader(){}

bool TachyonReader::open(void){
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

bool TachyonReader::nextBlock(){
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

TachyonReader::block_entry_type TachyonReader::getBlock(){
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

const int TachyonReader::has_format_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getFormatField(field_name, match)){
		U32 target = -1;
		for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
			//std::cerr << i << '/' << this->block.index_entry.n_format_streams << '\t' << this->block.index_entry.format_offsets[i].key << '\t' << this->header.entries[this->block.index_entry.format_offsets[i].key].ID << std::endl;
			if(this->block.index_entry.format_offsets[i].key == match->IDX){
				target = i;
				break;
			}
		}
		//std::cerr << "target stream is: " << target << std::endl;
		return(target);
	}
	return(-2);
}

const int TachyonReader::has_info_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getInfoField(field_name, match)){
		U32 target = -1;
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			//std::cerr << i << '/' << this->block.index_entry.n_info_streams << '\t' << this->block.index_entry.info_offsets[i].key << '\t' << this->header.entries[this->block.index_entry.info_offsets[i].key].ID << std::endl;
			if(this->block.index_entry.info_offsets[i].key == match->IDX){
				target = i;
				break;
			}
		}
		return(target);
	}
	return(-2);
}

const int TachyonReader::has_filter_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getFilterField(field_name, match)){
		U32 target = -1;
		for(U32 i = 0; i < this->block.index_entry.n_filter_streams; ++i){
			if(this->block.index_entry.filter_offsets[i].key == match->IDX){
				target = i;
				break;
			}
		}
		return(target);
	}
	return(-2);
}

const std::vector<bool> TachyonReader::get_info_field_pattern_matches(const std::string& field_name) const{
	int local_info_field_id = this->has_info_field(field_name);
	std::vector<bool> ret;
	if(local_info_field_id >= 0){
		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		ret.resize(this->block.index_entry.n_info_patterns, false);
		for(U32 i = 0; i < this->block.index_entry.n_info_patterns; ++i){
			//std::cerr << i << '\t' << this->block.index_entry.info_bit_vectors[i][local_info_field_id] << std::endl;
			ret[i] = this->block.index_entry.info_bit_vectors[i][local_info_field_id];
		}
	}
	return(ret);
}

const std::vector<bool> TachyonReader::get_format_field_pattern_matches(const std::string& field_name) const{
	int local_format_field_id = this->has_format_field(field_name);
	std::vector<bool> ret;
	if(local_format_field_id >= 0){
		// Collect all matches
		// Place in array
		// 0 = false, 1 = true
		ret.resize(this->block.index_entry.n_format_patterns, false);
		for(U32 i = 0; i < this->block.index_entry.n_format_patterns; ++i){
			//std::cerr << i << '\t' << this->block.index_entry.format_bit_vectors[i][local_format_field_id] << std::endl;
			ret[i] = this->block.index_entry.format_bit_vectors[i][local_format_field_id];
		}
	}
	return(ret);
}

}
