#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "TomahawkImportWriter.h"

namespace Tomahawk {

TomahawkImportWriter::TomahawkImportWriter() :
	flush_limit(1000000),
	n_variants_limit(1024),
	n_blocksWritten(0),
	n_variants_written(0),
	n_variants_complex_written(0),
	largest_uncompressed_block(0),
	buffer_encode_rle(flush_limit*2),
	buffer_encode_simple(flush_limit*2),
	buffer_meta(flush_limit*10), // meta joins all other buffers
	buffer_metaComplex(flush_limit*2),
	filter_hash_pattern(5012),
	info_hash_pattern(5012),
	info_hash_streams(5012),
	format_hash_pattern(5012),
	format_hash_streams(5012),
	filter_hash_streams(5012),
	buffer_ppa(100000),
	buffer_general(300000),
	encoder(nullptr),
	vcf_header(nullptr),
	info_containers(new stream_container[100]),
	format_containers(new stream_container[100])
{
	for(U32 i = 0; i < 100; ++i){
		this->info_containers[i].resize(65536*4);
		this->format_containers[i].resize(65536*1000);
	}
}

TomahawkImportWriter::~TomahawkImportWriter(){
	delete this->encoder;
	this->buffer_encode_rle.deleteAll();
	this->buffer_meta.deleteAll();
	this->buffer_encode_simple.deleteAll();
	this->buffer_metaComplex.deleteAll();
	this->buffer_general.deleteAll();
}

bool TomahawkImportWriter::Open(const std::string output){
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

void TomahawkImportWriter::WriteFinal(void){

}

void TomahawkImportWriter::setHeader(VCF::VCFHeader& header){
	this->vcf_header = &header;
	this->encoder = new encoder_type(header.samples);
}


bool TomahawkImportWriter::add(bcf_entry_type& entry, const U32* const ppa){
	// Keep positions
	// If the entry needs to be filtered out
	// then we roll back to these positions
	// In practice we simply move the pointer back
	const U64 meta_start_pos = this->buffer_meta.pointer;
	const U64 simple_start_pos = this->buffer_encode_simple.pointer;
	const U64 rle_start_pos  = this->buffer_encode_rle.pointer;
	meta_base_type meta;

	// Perform run-length encoding
	U64 n_runs = 0;
	std::cerr << "before encode" << std::endl;
	/*
	if(!this->encoder->Encode(entry, meta, this->buffer_encode_rle, this->buffer_encode_simple, n_runs, ppa)){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		return false;
	}
	*/
	std::cerr << "after encode" << std::endl;

	// Parse BCF
	this->parseBCF(meta, entry);
	std::cerr << "after parse bcf" << std::endl;

	// If the current minPosition is 0
	// then this is the first entry we've seen
	// in this contig. Keep the current position
	// as the last one we've seen
	if(this->totempole_entry.minPosition == 0)
		this->totempole_entry.minPosition = entry.body->POS + 1;

	// Update max position
	this->totempole_entry.maxPosition = entry.body->POS + 1;

	// Push meta to buffer
	// update complex offset position
	meta.virtual_offset_cold_meta = this->buffer_metaComplex.pointer;
	this->buffer_meta += meta;

	// RLE using this word size
	U32 w = ceil(ceil(log2(this->vcf_header->samples + 1))/8);
	if((w > 2) & (w < 4)) w = 4;
	else if(w > 4) w = 8;

	switch(w){
	case 1: this->buffer_meta += (BYTE)n_runs; break;
	case 2: this->buffer_meta += (U16)n_runs; break;
	case 4: this->buffer_meta += (U32)n_runs; break;
	case 8: this->buffer_meta += (U64)n_runs; break;
	default:
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word-size!" << std::endl;
		exit(1); // unrecoverable error
	}

	// Complex meta data
	Core::EntryColdMeta test;
	//if(!test.write(entry, this->buffer_metaComplex)){
	//	std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
	//	return false;
	//}

	// Update number of entries in block
	++this->totempole_entry.n_variants;

	return true;
}

bool TomahawkImportWriter::parseBCF(meta_base_type& meta, bcf_entry_type& entry){
	//std::cerr << this->body->CHROM << ':' << this->body->POS+1 << '\t' << this->body->n_allele << '\t' << this->body->n_fmt << '\t' << this->body->n_info << '\t' << this->body->n_sample << std::endl;
	U32 internal_pos = entry.filter_start;

	// At FILTER
	// Typed vector
	const bcf_entry_type::base_type& filter_key = *reinterpret_cast<const bcf_entry_type::base_type* const>(&entry.data[internal_pos++]);
	U32 n_filter = filter_key.high;
	if(n_filter == 15) n_filter = entry.getInteger(filter_key.low, internal_pos);
	entry.n_filter = n_filter;
	entry.filter_key = filter_key;

	S32 val = 0;
	while(entry.nextFilter(val, internal_pos)){
		// Hash FILTER value
		// Filter fields have no values
		U32* hash_map_ret = nullptr;
		U32 temp = val;
		if(this->filter_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			// exists
		} else {
			U32 tot = this->filter_values.size();
			this->filter_hash_streams.SetItem(&temp, tot, sizeof(U32));
			this->filter_values.push_back(val);
		}
	}

	// At INFO
	U32 info_length;
	BYTE info_value_type;
	while(entry.nextInfo(val, info_length, info_value_type, internal_pos)){
		// Hash INFO values
		U32* hash_map_ret = nullptr;
		U32 mapID = 0;
		U32 temp = val;
		if(this->info_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			mapID = *hash_map_ret;
		} else {
			U32 tot = this->info_values.size();
			this->info_hash_streams.SetItem(&temp, tot, sizeof(U32));
			mapID = tot;
			this->info_values.push_back(val);
		}

		//
		stream_container& target_container = this->info_containers[mapID];
		if(this->info_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(4);
			else if(info_value_type == 1) target_container.setType(4);
			else if(info_value_type == 2) target_container.setType(4);
			else if(info_value_type == 3) target_container.setType(4);
			else if(info_value_type == 5) target_container.setType(7);
			else if(info_value_type == 7) target_container.setType(0);
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getInteger(info_value_type, internal_pos);
			}
		}
		// Floats
		else if(info_value_type == 5){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_value_type == 7){
			for(U32 j = 0; j < info_length; ++j){
				target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)info_value_type << std::endl;
			exit(1);
		}
	}

#if BCF_ASSERT == 1
	// Assert all FILTER and INFO data have been successfully
	// parsed. This is true when the byte pointer equals the
	// start position of the FORMAT fields which are encoded
	// in the meta header structure
	assert(internal_pos == (entry.body->l_shared + sizeof(U32)*2));
#endif

	while(entry.nextFormat(val, info_length, info_value_type, internal_pos)){
		// Hash FORMAT values
		U32* hash_map_ret = nullptr;
		U32 mapID = 0;
		U32 temp = val;
		if(this->format_hash_streams.GetItem(&temp, hash_map_ret, sizeof(U32))){
			mapID = *hash_map_ret;
		} else {
			U32 tot = this->format_values.size();
			this->format_hash_streams.SetItem(&temp, tot, sizeof(U32));
			mapID = tot;
			this->format_values.push_back(val);
			std::cerr << Helpers::timestamp("DEBUG") << val << '\t' << info_length << '\t' << (U32)info_value_type << std::endl;
		}
		std::cerr << "here" << std::endl;

		std::cerr << "mapid: " << mapID << std::endl;
		std::cerr << this->format_containers << std::endl;
		std::cerr << &this->format_containers[1] << std::endl;
		std::cerr << this->format_containers[0].n_entries << std::endl;
		if(this->format_containers[mapID].n_entries == 0){
			this->format_containers[mapID].setStrideSize(info_length);
		}
		std::cerr << "after update" << std::endl;

		if(mapID == 0){
			switch(info_value_type){
			case 1: internal_pos += this->vcf_header->samples * sizeof(SBYTE) * info_length; break;
			case 2: internal_pos += this->vcf_header->samples * sizeof(S16)   * info_length; break;
			case 3: internal_pos += this->vcf_header->samples * sizeof(S32)   * info_length; break;
			}
			continue;
		}
		std::cerr << "after update2" << std::endl;

		stream_container& target_container = this->format_containers[mapID];
		if(this->format_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(4);
			else if(info_value_type == 1) target_container.setType(4);
			else if(info_value_type == 2) target_container.setType(4);
			else if(info_value_type == 3) target_container.setType(4);
			else if(info_value_type == 5) target_container.setType(7);
			else if(info_value_type == 7) target_container.setType(0);
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getInteger(info_value_type, internal_pos);
			}
		}
		// Floats
		else if(info_value_type == 5){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_value_type == 7){
			for(U32 j = 0; j < this->vcf_header->samples*info_length; ++j){
				target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)info_value_type << std::endl;
			exit(1);
		}
	}

	// Hash FILTER pattern
	U32 mapID = 0;
	U32* hash_map_ret = nullptr;
	const U64 hash_filter_vector = entry.hashFilter();
	if(this->filter_hash_pattern.GetItem(&hash_filter_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		U32 tot = this->filter_patterns.size();
		this->filter_hash_pattern.SetItem(&hash_filter_vector, tot, sizeof(U64));
		this->filter_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.filterPointer; ++i){
			this->filter_patterns[tot].push_back(entry.filterID[i]);
		}
		assert(tot < 65536);
		mapID = tot;
	}
	// Store this map in the meta
	meta.FILTER_map_ID = mapID;

	// Hash INFO pattern
	mapID = 0;
	hash_map_ret = nullptr;
	// Hash INFO vector of identifiers
	const U64 hash_info_vector = entry.hashInfo();
	if(this->info_hash_pattern.GetItem(&hash_info_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		U32 tot = this->info_patterns.size();
		this->info_hash_pattern.SetItem(&hash_info_vector, tot, sizeof(U64));
		this->info_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.infoPointer; ++i){
			this->info_patterns[tot].push_back(entry.infoID[i]);
		}
		assert(tot < 65536);
		mapID = tot;
	}
	// Store this map in the meta
	meta.INFO_map_ID = mapID;

	// Hash FORMAT pattern
	mapID = 0;
	hash_map_ret = nullptr;
	const U64 hash_format_vector = entry.hashFormat();
	if(this->format_hash_pattern.GetItem(&hash_format_vector, hash_map_ret, sizeof(U64))){
		mapID = *hash_map_ret;
	} else {
		U32 tot = this->format_patterns.size();
		this->format_hash_pattern.SetItem(&hash_format_vector, tot, sizeof(U64));
		format_patterns.push_back(std::vector<U32>());
		for(U32 i = 0; i < entry.formatPointer; ++i){
			this->format_patterns[tot].push_back(entry.formatID[i]);
		}
		assert(tot < 65536);
		mapID = tot;
	}
	// Store this map in the meta
	meta.FORMAT_map_ID = mapID;

	// Return
	return true;
}

// flush and write
bool TomahawkImportWriter::flush(const U32* ppa){
	if(this->buffer_meta.size() == 0)
		return false;

	// Update sizes of streams
	this->totempole_entry.byte_offset = this->streamTomahawk.tellp(); // IO offset in Tomahawk output
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

S32 TomahawkImportWriter::recodeStreamStride(stream_container& stream){
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
	this->buffer_general.reset();

	return(ret_size);
}

bool TomahawkImportWriter::checkUniformity(stream_container& stream){
	const U32 stride_size = stream.header.stride;
	if(stride_size == -1)
		return false;

	U32 stride_update = stride_size;

	switch(stream.header.controller.type){
	case 4: stride_update *= sizeof(S32);   break;
	case 7: stride_update *= sizeof(float); break;
	case 0: stride_update *= sizeof(char);  break;
	default: return false; break;
	}

	bool is_uniform = true;
	const U64 first_hash = XXH64(&stream.buffer_data.data[0], stride_update, 1337);
	for(U32 i = stride_update; i < stream.n_entries; i+= stride_update){
		if(XXH64(&stream.buffer_data.data[i], stride_update, 1337) != first_hash){
			std::cerr << Helpers::timestamp("DEBUG") << "Not uniform" << std::endl;
			is_uniform = false;
			break;
		}
	}
	std::cerr << Helpers::timestamp("DEBUG") << "Uniformity: " << is_uniform << std::endl;
	return(is_uniform);
}

S32 TomahawkImportWriter::recodeStream(stream_container& stream){
	S32 ret_size = -1;

	if(this->checkUniformity(stream)){
		std::cerr << "do stuff with me" << std::endl;
		stream.header.controller.uniform = true;
		stream.header.controller.mixedStride = false;
		stream.header.controller.encoder = 0;
	}

	if(stream.header.controller.type == 4){
		const S32* dat = reinterpret_cast<const S32*>(stream.buffer_data.data);
		S32 min = dat[0];
		S32 max = dat[0];
		S32 prev_value = dat[0];
		bool is_uniform = true;

		for(U32 j = 1; j < stream.n_entries; ++j){
			if(dat[j] < min) min = dat[j];
			if(dat[j] > max) max = dat[j];
			if(prev_value != dat[j]) is_uniform = false;
			prev_value = dat[j];
		}
		std::cerr << Helpers::timestamp("DEBUG") << "Uniformity: " << is_uniform << std::endl;

		BYTE byte_width = 0;
		if(min < 0) byte_width = ceil((ceil(log2(abs(min) + 1))+1)/8);  // One bit is used for sign
		else byte_width = ceil(ceil(log2(max + 1))/8);

		if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
		else if(byte_width > 4) byte_width = 8;
		if(byte_width == 0) byte_width = 1;

		// Phase 2
		// Here we re-encode values using the smallest possible
		// word-size
		// Todo: uniform over stride size
		if(is_uniform){
			// Non-negative
			if(min >= 0){
				switch(byte_width){
				case 1: this->buffer_general += (BYTE)min; break;
				case 2: this->buffer_general += (U16)min; break;
				case 4: this->buffer_general += (U32)min; break;
				case 8: this->buffer_general += (U64)min; break;
				default: std::cerr << "illegal: " << std::endl; exit(1);
				}
			} else {
				switch(byte_width){
				case 1: this->buffer_general += (SBYTE)min; break;
				case 2: this->buffer_general += (S16)min; break;
				case 4: this->buffer_general += (S32)min; break;
				default: std::cerr << "illegal" << std::endl; exit(1);
				}
			}

			std::cerr << "Return is uniform" << std::endl;
			// Write out data as a literal
			this->streamTomahawk << stream.header;
			this->streamTomahawk.write(this->buffer_general.data, this->buffer_general.pointer);
			return(this->buffer_general.pointer);
		} else {
			//std::cerr << "non-uniform" << std::endl;
			dat = reinterpret_cast<const S32*>(stream.buffer_data.data);
			// Is non-negative
			if(min >= 0){
				if(byte_width == 1){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (BYTE)*(dat++);
				} else if(byte_width == 2){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U16)*(dat++);
				} else if(byte_width == 4){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U32)*(dat++);
				} else if(byte_width == 8){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (U64)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
			// Is negative
			else {
				if(byte_width == 1){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (SBYTE)*(dat++);
				} else if(byte_width == 2){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (S16)*(dat++);
				} else if(byte_width == 4){
					for(U32 j = 0; j < stream.n_entries; ++j)
						this->buffer_general += (S32)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
		}

		this->gzip_controller.Deflate(this->buffer_general);
		this->streamTomahawk << stream.header;
		this->streamTomahawk << this->gzip_controller;
		ret_size = this->gzip_controller.buffer.size();
		this->gzip_controller.Clear();
		this->buffer_general.reset();
		return(ret_size);
	}
	// Is not an integer
	else {
		this->gzip_controller.Deflate(stream.buffer_data);
		this->streamTomahawk << stream.header;
		this->streamTomahawk << this->gzip_controller;
		ret_size = this->gzip_controller.buffer.size();
		this->gzip_controller.Clear();
	}

	return(ret_size);
}

void TomahawkImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &Constants::OUTPUT_SUFFIX[0], Constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tomahawk */
