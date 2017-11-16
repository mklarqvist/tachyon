#include <fstream>
#include "TomahawkImporter.h"

namespace Tomahawk {

TomahawkImporter::TomahawkImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint) :
	checkpoint_size(checkpoint),
	block_flush_limit(65536),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader_(inputFile),
	writer_(),
	header_(nullptr)
{}

TomahawkImporter::~TomahawkImporter(){}

bool TomahawkImporter::Build(){
	std::ifstream temp(this->inputFile, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT")  << "Failed to open file (" << this->inputFile << ")..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if((BYTE)tempData[0] == IO::Constants::GZIP_ID1 && (BYTE)tempData[1] == IO::Constants::GZIP_ID2){
		if(!this->BuildBCF()){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool TomahawkImporter::BuildBCF(void){
	bcf_reader_type reader;
	if(!reader.open(this->inputFile)){
		std::cerr << Helpers::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	this->header_ = &reader.header;
	if(this->header_->samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "No samples detected in header..." << std::endl;
		return false;
	}

	if(this->header_->samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller and update PPA controller
	this->encoder.setSamples(this->header_->samples);
	this->permutator.setSamples(this->header_->samples);

	this->writer_.setHeader(reader.header);
	if(!this->writer_.Open(this->outputPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	// Resize containers
	const U32 resize_to = this->checkpoint_size * sizeof(U32) * this->header_->samples;
	this->meta_hot_container.resize(resize_to);
	this->meta_cold_container.resize(resize_to);
	this->gt_rle_container.resize(resize_to);
	this->gt_simple_container.resize(resize_to);

	while(true){
		if(!reader.getVariants(this->checkpoint_size))
			break;

		// Reset permutate
		if(!this->permutator.build(reader)){
			std::cerr << "fail" << std::endl;
			return false;
		}

		std::cerr << "After permuter" << std::endl;
		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->parseBCFLine(reader[i])){
				std::cerr << "failed to parse" << std::endl;
				return false;
			}
		}

		std::cerr << "After parse" << std::endl;
		++this->header_->getContig(reader[0].body->CHROM); // update block count for this contigID
		this->writer_.flush(this->permutator.getPPA());
		this->writer_.n_variants_written += reader.size();

		// Reset permutator
		this->permutator.reset();
	}

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_->getContig(*this->sort_order_helper.contigID);
	this->writer_.flush(this->permutator.getPPA());
	this->writer_.WriteFinal();
	this->writer_.n_variants_written += reader.size();

	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;
	return(true);
}

bool TomahawkImporter::parseBCFLine(bcf_entry_type& entry){
	// Assert position is in range
	if(entry.body->POS + 1 > this->header_->getContig(entry.body->CHROM).length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << (*this->header_)[entry.body->CHROM].name << ':' << entry.body->POS+1 << " > reported max size of contig (" << (*this->header_)[entry.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Perform run-length encoding
	U64 n_runs = 0;
	std::cerr << "before encode" << std::endl;

	meta_type meta;
	if(!this->encoder.Encode(entry, meta, this->gt_rle_container, this->gt_simple_container, n_runs, this->permutator.ppa))
		return false;

	U32 test = 24;
	U32 ret = this->info_fields.setGet(test);
	std::cerr << "got stuff back: "  << ret << std::endl;

	/*
	// Keep positions
	// If the entry needs to be filtered out
	// then we roll back to these positions
	// In practice we simply move the pointer back
	const U64 meta_start_pos   = this->buffer_meta.pointer;
	const U64 simple_start_pos = this->buffer_encode_simple.pointer;
	const U64 rle_start_pos    = this->buffer_encode_rle.pointer;
	meta_base_type meta;

	// Perform run-length encoding
	U64 n_runs = 0;
	std::cerr << "before encode" << std::endl;
	if(!this->encoder->Encode(entry, meta, this->buffer_encode_rle, this->buffer_encode_simple, n_runs, ppa)){
		this->buffer_meta.pointer = meta_start_pos; // roll back
		this->buffer_encode_rle.pointer  = rle_start_pos; // roll back
		this->buffer_encode_simple.pointer = simple_start_pos;
		return false;
	}
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
	Support::EntryColdMeta test;
	if(!test.write(entry, this->buffer_metaComplex)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	// Update number of entries in block
	++this->totempole_entry.n_variants;
	*/

	return true;
}

bool TomahawkImporter::parseBCFBody(meta_type& meta, bcf_entry_type& entry){
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
		this->filter_fields.set(val);
	}

	// At INFO
	U32 info_length;
	BYTE info_value_type;
	while(entry.nextInfo(val, info_length, info_value_type, internal_pos)){
		// Hash INFO values
		const U32 mapID = this->info_fields.setGet(val);

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

} /* namespace Tomahawk */
