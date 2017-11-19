#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "ImportWriter.h"

namespace Tachyon {

ImportWriter::ImportWriter() :
	n_blocksWritten(0),
	n_variants_written(0)
{
}

ImportWriter::~ImportWriter(){}

bool ImportWriter::Open(const std::string output){
	this->filename = output;
	this->CheckOutputNames(output);
	this->streamTomahawk.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->streamTomahawk.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

void ImportWriter::WriteFinal(void){

}

// flush and write
bool ImportWriter::flush(void){
	//if(this->buffer_meta.size() == 0)
	//	return false;

	// Update sizes of streams
	//this->totempole_entry.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	//this->totempole_entry.l_meta = this->buffer_meta.pointer;
	//this->totempole_entry.l_meta_complex = this->buffer_metaComplex.pointer;
	//this->totempole_entry.l_gt_rle = this->buffer_encode_rle.pointer;
	//this->totempole_entry.l_gt_simple = this->buffer_encode_simple.pointer;

	/*
	Index::IndexBlockEntry block_entry;
	std::cerr << "ID: " << block_entry.contigID << std::endl;
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO: " << this->info_values.size() << " entries -> " << ceil((float)this->info_values.size()/8) << " bytes" << std::endl;
	std::cerr << "External before: " << block_entry.info_bit_vectors << std::endl;
	std::cerr << &block_entry << std::endl;
	block_entry.constructBitVector(Index::IndexBlockEntry::INDEX_BLOCK_TARGET::INDEX_INFO, this->info_hash_streams, this->info_values, this->info_patterns);
	std::cerr << "External after: " << block_entry.info_bit_vectors << std::endl;

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT: " << this->format_values.size() << " entries -> " << ceil((float)this->format_values.size()/8) << " bytes" << std::endl;
	block_entry.constructBitVector(Index::IndexBlockEntry::INDEX_BLOCK_TARGET::INDEX_FORMAT, this->format_hash_streams, this->format_values, this->format_patterns);

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FILTER: " << this->filter_values.size() << " entries -> " << ceil((float)this->filter_values.size()/8) << " bytes" << std::endl;
	block_entry.constructBitVector(Index::IndexBlockEntry::INDEX_BLOCK_TARGET::INDEX_FILTER, this->filter_hash_streams, this->filter_values, this->filter_patterns);


	std::cerr << block_entry.info_bit_vectors[0][6] << std::endl;
*/

	/*
	// Split U32 values into 4 streams
	const U32 partition = this->vcf_header->samples;
	//buffer_type test(partition*sizeof(U32)*10);
	bytePreprocessor(ppa, partition, this->buffer_general.data);

	// test bitshuffle
	memset(this->buffer_ppa.data, 0, 2*sizeof(BYTE)*partition);
	bytePreprocessBits(&this->buffer_general.data[partition*2], partition, this->buffer_ppa.data);
	bytePreprocessBits(&this->buffer_general.data[partition*3], partition, &this->buffer_ppa.data[partition]);
	this->buffer_ppa.pointer = 2*partition;
	this->gzip_controller.Deflate(this->buffer_ppa);
	this->streamTomahawk << this->gzip_controller;
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "PPA\t" << 2*partition << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)2*partition/this->gzip_controller.buffer.size() << std::endl;
	//std::cout << "PPA\t0\t" << partition*2 << '\t' << this->gzip_controller.buffer.size() << '\n';
	this->gzip_controller.Clear();
	this->buffer_general.reset();

	// Merge data to single buffer
	// and compress
	//
	// Meta
	//this->buffer_meta += this->buffer_metaComplex;
	this->gzip_controller.Deflate(this->buffer_meta); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "META\t" << this->buffer_meta.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_meta.size() / this->gzip_controller.buffer.size() << std::endl;
	//std::cout << "META\t0\t" << this->buffer_meta.size() << '\t' << this->gzip_controller.buffer.size() << std::endl;
	this->gzip_controller.Clear(); // Clean up gzip controller

	this->gzip_controller.Deflate(this->buffer_metaComplex); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "META-C\t" << this->buffer_metaComplex.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_metaComplex.size() / this->gzip_controller.buffer.size() << std::endl;
	//std::cout << "META-C\t0\t" << this->buffer_metaComplex.size() << '\t' << this->gzip_controller.buffer.size() << std::endl;
	this->gzip_controller.Clear(); // Clean up gzip controller

	// RLE
	//this->buffer_encode_rle += this->buffer_encode_simple;
	this->gzip_controller.Deflate(this->buffer_encode_rle); // Deflate block
	this->streamTomahawk << this->gzip_controller; // Write tomahawk output
	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "RLE\t" << this->buffer_encode_rle.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_encode_rle.size() / this->gzip_controller.buffer.size() << std::endl;
	//std::cout << "RLE\t0\t" << this->buffer_encode_rle.size() << '\t' << this->gzip_controller.buffer.size() << std::endl;
	this->gzip_controller.Clear(); // Clean up gzip controller

	// If there is any simple data
	if(this->buffer_encode_simple.size() > 0){
		this->gzip_controller.Deflate(this->buffer_encode_simple); // Deflate block
		this->streamTomahawk << this->gzip_controller; // Write tomahawk output
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "RLE-S\t" << this->buffer_encode_simple.size() << '\t' << this->gzip_controller.buffer.size() << '\t' << (float)this->buffer_encode_simple.size() / this->gzip_controller.buffer.size() << std::endl;
		//std::cout << "RLE-S\t0\t" << this->buffer_encode_simple.size() << '\t' << this->gzip_controller.buffer.size() << std::endl;
		this->gzip_controller.Clear(); // Clean up gzip controller
	}

	// Keep track of largest block observed
	if(this->buffer_encode_rle.size() > this->largest_uncompressed_block)
		this->largest_uncompressed_block = this->buffer_encode_rle.size();

	//this->totempole_entry.l_uncompressed = this->buffer_meta.size(); // Store uncompressed size
	this->totempole_entry.byte_offset_end = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
	++this->n_blocksWritten; // update number of blocks written
	this->n_variants_written += this->totempole_entry.n_variants; // update number of variants written
	this->totempole_entry.reset();

	// Dispatch values into streams
	for(U32 i = 0; i < this->info_values.size(); ++i){
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "Recoding: INFO " << this->info_values[i] << " with stride: " << this->info_containers[i].header.stride << std::endl;
		S32 ret_size = this->recodeStream(this->info_containers[i]);
		if(ret_size == -1){
			std::cerr << "failed recode @ " << i << std::endl;
			exit(1);
		}
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO " << this->info_values[i] << '\t' << this->info_containers[i].buffer_data.size() << '\t' << ret_size << '\t' << (float)this->info_containers[i].buffer_data.size()/ret_size << '\t' << this->info_containers[i].header.stride << std::endl;
		//std::cout << "INFO\t" << this->info_values[i] << '\t' << this->info_containers[i].buffer_data.size() << '\t' << ret_size << std::endl;
		if(this->info_containers[i].header.stride == -1){
			S32 ret_stride_size = this->recodeStreamStride(this->info_containers[i]);
			std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "INFO-ADD " << this->info_values[i] << '\t' << this->info_containers[i].buffer_strides.size() << '\t' << ret_stride_size << '\t' << (float)this->info_containers[i].buffer_strides.size()/ret_stride_size << std::endl;
			//std::cout << "INFO-ADD\t" << this->info_values[i] << '\t' << this->info_containers[i].buffer_strides.size() << '\t' << ret_stride_size << std::endl;
		}
	}

	// Dispatch values into streams
	for(U32 i = 0; i < this->format_values.size(); ++i){
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "Recoding: FORMAT " << this->format_values[i] << " with stride: " << this->format_containers[i].header.stride << std::endl;
		S32 ret_size = this->recodeStream(this->format_containers[i]);
		if(ret_size == -1){
			std::cerr << "failed recode @ " << i << std::endl;
			exit(1);
		}
		std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT " << this->format_values[i] << '\t' << this->format_containers[i].buffer_data.size() << '\t' << ret_size << '\t' << (float)this->format_containers[i].buffer_data.size()/ret_size << '\t' << this->format_containers[i].header.stride << std::endl;
		//std::cout << "FORMAT\t" << this->format_values[i] << '\t' << this->format_containers[i].buffer_data.size() << '\t' << ret_size << std::endl;
		if(this->format_containers[i].header.stride == -1){
			S32 ret_stride_size = this->recodeStreamStride(this->format_containers[i]);
			std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "FORMAT-ADD " << this->format_values[i] << '\t' << this->format_containers[i].buffer_strides.size() << '\t' << ret_stride_size << '\t' << (float)this->format_containers[i].buffer_strides.size()/ret_stride_size << std::endl;
			//std::cout << "FORMAT-ADD\t" << this->format_values[i] << '\t' << this->format_containers[i].buffer_strides.size() << '\t' << ret_stride_size << std::endl;
		}
	}

	this->reset(); // reset buffers

	std::cerr << Helpers::timestamp("DEBUG","FLUSH") << "END FLUSH" << std::endl;

	return true;
}

S32 ImportWriter::recodeStreamStride(stream_container& stream){
	S32 ret_size = -1;
	const U32* const strides = reinterpret_cast<const U32* const>(stream.buffer_strides.data);
	if(stream.buffer_strides.size() % sizeof(U32) != 0){
		std::cerr << Helpers::timestamp("ERROR") << "Buffer is truncated!" << std::endl;
		return -1;
	}
	const U32 n_entries = stream.buffer_strides.size() / sizeof(U32);
	U32 min = strides[0];
	U32 max = strides[0];
	U32 prev_value = strides[0];
	bool is_uniform = true;

	for(U32 i = 0; i < n_entries; ++i){
		if(strides[i] < min) min = strides[i];
		if(strides[i] > max) max = strides[i];
		if(prev_value != strides[i]) is_uniform = false;
		prev_value = strides[i];
	}

	BYTE byte_width = byte_width = ceil(ceil(log2(max + 1)) / 8 );
	if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if(byte_width > 4){
		std::cerr << Helpers::timestamp("ERROR") << "Too large" << std::endl;
		return(-1);
	}

	if(byte_width == 1){
		for(U32 j = 0; j < n_entries; ++j)
			this->buffer_general += (BYTE)strides[j];
	} else if(byte_width == 2){
		for(U32 j = 0; j < n_entries; ++j)
			this->buffer_general += (U16)strides[j];
	} else if(byte_width == 4){
		for(U32 j = 0; j < n_entries; ++j)
			this->buffer_general += (U32)strides[j];
	} else {
		std::cerr << "illegal" << std::endl;
		exit(1);
	}

	this->gzip_controller.Deflate(this->buffer_general);
	this->streamTomahawk << stream.header_stride;
	this->streamTomahawk << this->gzip_controller;
	ret_size = this->gzip_controller.buffer.size();
	this->gzip_controller.Clear();
	this->buffer_gener
	*al.reset();

	return(ret_size);
	*/
}


void ImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &Constants::OUTPUT_SUFFIX[0], Constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tomahawk */
