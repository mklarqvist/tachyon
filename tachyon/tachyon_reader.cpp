#include "tachyon_reader.h"

namespace tachyon{

bool TachyonReader::open(void){
	if(this->input_file.size() == 0){
		std::cerr << "no filename" << std::endl;
		return false;
	}

	this->stream.open(this->input_file, std::ios::binary | std::ios::in | std::ios::ate);
	this->filesize = (U64)this->stream.tellg();

	if(!this->stream.good()){
		std::cerr << "failed to read file" << std::endl;
		return false;
	}

	this->stream.seekg(0);
	if(!this->stream.good()){
		std::cerr << "failed to rewind" << std::endl;
		return false;
	}

	// Load header
	this->stream << this->header;

	if(!this->stream.good()){
		std::cerr << "failed to get header" << std::endl;
		return false;
	}

	// Keep track of start position
	const U64 return_pos = this->stream.tellg();

	// Seek to EOF and make check
	this->stream.seekg(this->filesize - 32);
	BYTE eof_data[32];
	utility::HexToBytes(constants::TACHYON_FILE_EOF, &eof_data[0]);
	BYTE eof_match[32];
	this->stream.read((char*)&eof_match[0], 32);
	for(U32 i = 0; i < 32; ++i){
		if(eof_data[i] != eof_match[i]){
			std::cerr << "File is truncated!" << std::endl;
			return false;
		}
	}

	// Seek back to find start of index
	// Seek to that start of index
	// Load index
	// Seek back to start of the file
	this->stream.seekg(this->filesize - 32 - sizeof(U64));
	this->stream.read((char*)reinterpret_cast<char*>(&this->l_data), sizeof(U64));
	this->stream.seekg(this->l_data);
	this->stream >> this->index;
	this->stream.seekg(return_pos);

	return(this->stream.good());
}

bool TachyonReader::get_next_block(){
	// If the stream is faulty then return
	if(!this->stream.good()){
		std::cerr << "faulty stream" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->stream.tellg() == this->l_data)
		return false;

	// Reset and re-use
	this->block.clear();

	// Attempts to read a YON block with the provided
	// settings
	if(!this->block.read(stream, this->settings))
		return false;

	// Internally decompress available data
	if(!this->codec_manager.decompress(this->block))
		return false;

	// All passed
	return true;
}

const int TachyonReader::has_format_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getEntry(field_name, match)){
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
	if(this->header.getEntry(field_name, match)){
		U32 target = -1;
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			//std::cerr << i << '/' << this->block.index_entry.n_info_streams << '\t' << this->block.index_entry.info_offsets[i].key << '\t' << this->header.entries[this->block.index_entry.info_offsets[i].key].ID << std::endl;
			if(this->block.index_entry.info_offsets[i].key == match->IDX){
				target = i;
				break;
			}
		}
		//std::cerr << "target stream is: " << target << std::endl;
		return(target);
	}
	return(-2);
}

const int TachyonReader::has_filter_field(const std::string& field_name) const{
	core::HeaderMapEntry* match = nullptr;
	if(this->header.getEntry(field_name, match)){
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
