#include <fstream>
#include "Importer.h"

namespace Tomahawk {

Importer::Importer(std::string inputFile, std::string outputPrefix, const U32 checkpoint) :
	checkpoint_size(checkpoint),
	block_flush_limit(65536),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader_(inputFile),
	writer_(),
	header_(nullptr),
	info_containers(new stream_container[100]),
	format_containers(new stream_container[100]),
	filter_containers(new stream_container[100])
{
	for(U32 i = 0; i < 100; ++i){
		this->info_containers[i].resize(65536*4);
		this->format_containers[i].resize(65536*4);
		this->filter_containers[i].resize(65536*4);
	}
}

Importer::~Importer(){}

void Importer::resetHashes(void){
	this->info_fields.clear();
	this->info_patterns.clear();
	this->format_fields.clear();
	this->format_patterns.clear();
	this->filter_fields.clear();
	this->filter_patterns.clear();
}

void Importer::resetContainers(void){
	for(U32 i = 0; i < 100; ++i){
		this->info_containers[i].reset();
		this->format_containers[i].reset();
		this->filter_containers[i].reset();
	}
	this->meta_hot_container.reset();
	this->meta_cold_container.reset();
	this->gt_rle_container.reset();
	this->gt_simple_container.reset();
}

bool Importer::Build(){
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

bool Importer::BuildBCF(void){
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

	IO::TGZFController cont;
	cont.buffer.resize(256000);

	while(true){
		if(!reader.getVariants(this->checkpoint_size))
			break;

		// Reset permutate
		std::cerr << "reader: " << reader.size() << std::endl;
		if(!this->permutator.build(reader)){
			std::cerr << "fail permutator" << std::endl;
			return false;
		}

		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->parseBCFLine(reader[i])){
				std::cerr << "failed to parse" << std::endl;
				return false;
			}
		}

		//++this->header_->getContig(reader[0].body->CHROM); // update block count for this contigID
		//this->writer_.flush();
		//this->writer_.n_variants_written += reader.size();
		cont.Deflate(this->meta_hot_container.buffer_data);
		this->writer_.streamTomahawk << cont;
		std::cerr <<"META-HOT\t" << this->meta_hot_container.buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
		cont.Clear();

		cont.Deflate(this->meta_cold_container.buffer_data);
		this->writer_.streamTomahawk << cont;
		std::cerr <<"META-COLD\t" << this->meta_cold_container.buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
		cont.Clear();

		cont.Deflate(this->gt_rle_container.buffer_data);
		this->writer_.streamTomahawk << cont;
		std::cerr <<"GT RLE\t" << this->gt_rle_container.buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
		cont.Clear();

		cont.Deflate(this->gt_simple_container.buffer_data);
		this->writer_.streamTomahawk << cont;
		std::cerr <<"GT SIMPLE\t" << this->gt_simple_container.buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
		cont.Clear();


		// Reset permutator

		std::cerr << "PATTERNS: " << this->info_patterns.size() << '\t' << this->format_patterns.size() << '\t' << this->filter_patterns.size() << std::endl;
		std::cerr << "VALUES: " << this->info_fields.size() << '\t' << this->format_fields.size() << '\t' << this->filter_fields.size() << std::endl;
		std::cerr << "INFO: " << std::endl;
		for(U32 i = 0; i < this->info_fields.size(); ++i){
			if(this->info_containers[i].buffer_data.size() == 0)
				continue;

			this->info_containers[i].checkUniformity();
			this->info_containers[i].reformatStream();

			cont.Deflate(this->info_containers[i].buffer_data);
			this->writer_.streamTomahawk << cont;
			std::cerr << this->info_fields[i] << '\t' << this->info_containers[i].buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
			cont.Clear();

			if(this->info_containers[i].header.stride != -1){
				cont.Deflate(this->info_containers[i].buffer_strides);
				this->writer_.streamTomahawk << cont;
				std::cerr << this->info_fields[i] << "-ADD\t" << this->info_containers[i].buffer_strides.size() << '\t' << cont.buffer.size() << std::endl;
				cont.Clear();
			}
		}

		std::cerr << "FORMAT: " << std::endl;
		for(U32 i = 0; i < this->format_fields.size(); ++i){
			if(this->format_containers[i].buffer_data.size() == 0)
				continue;

			this->format_containers[i].checkUniformity();

			cont.Deflate(this->format_containers[i].buffer_data);
			this->writer_.streamTomahawk << cont;
			std::cerr << this->format_fields[i] << '\t' << this->format_containers[i].buffer_data.size() << '\t' << cont.buffer.size() << std::endl;
			cont.Clear();
			if(this->format_containers[i].header.stride != -1){
				cont.Deflate(this->format_containers[i].buffer_strides);
				this->writer_.streamTomahawk << cont;
				std::cerr << this->format_fields[i] << "-ADD\t" << this->format_containers[i].buffer_strides.size() << std::endl;
				cont.Clear();
			}
		}

		//
		//Index::IndexBlockEntry i;
		//i.constructBitVector(Index::IndexBlockEntry::INDEX_INFO, this->info_fields, this->info_patterns);
		//i.constructBitVector(Index::IndexBlockEntry::INDEX_FORMAT, this->format_fields, this->format_patterns);
		//i.constructBitVector(Index::IndexBlockEntry::INDEX_FILTER, this->filter_fields, this->filter_patterns);

		std::cerr << "clear before" << std::endl;
		this->permutator.reset();
		this->resetHashes();
		this->resetContainers();
		std::cerr << "clear after" << std::endl;
	}
	std::cerr << "done main" << std::endl;

	/*
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
	 */
	std::cerr << "all done" << std::endl;
	return(true);
}

bool Importer::parseBCFLine(bcf_entry_type& entry){
	// Assert position is in range
	if(entry.body->POS + 1 > this->header_->getContig(entry.body->CHROM).length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << (*this->header_)[entry.body->CHROM].name << ':' << entry.body->POS+1 << " > reported max size of contig (" << (*this->header_)[entry.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Perform run-length encoding
	U64 n_runs = 0;

	meta_type meta;
	if(!this->encoder.Encode(entry, meta, this->gt_rle_container, this->gt_simple_container, n_runs, this->permutator.manager.get()))
		return false;

	if(!this->parseBCFBody(meta, entry)){
		std::cerr << "bad" << std::endl;
		return false;
	}

	// RLE using this word size
	U32 w = ceil(ceil(log2(this->header_->samples + 1))/8);
	if((w > 2) & (w < 4)) w = 4;
	else if(w > 4) w = 8;

	switch(w){
	case 1: this->meta_hot_container += (BYTE)n_runs; break;
	case 2: this->meta_hot_container += (U16)n_runs; break;
	case 4: this->meta_hot_container += (U32)n_runs; break;
	case 8: this->meta_hot_container += (U64)n_runs; break;
	default:
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word-size!" << std::endl;
		exit(1); // unrecoverable error
	}

	// Complex meta data
	Core::EntryColdMeta test;
	if(!test.write(entry, this->meta_cold_container)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	// Update number of entries in block
	++this->totempole_entry.n_variants;


	return true;
}

bool Importer::parseBCFBody(meta_type& meta, bcf_entry_type& entry){
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
			if(info_value_type == 0)      target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 1) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 2) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 3) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 5) target_container.setType(Core::CORE_TYPE::TYPE_FLOAT);
			else if(info_value_type == 7) target_container.setType(Core::CORE_TYPE::TYPE_8B);
			if(info_value_type != 5) target_container.header.controller.signedness = 1;
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		//std::cerr << val << '\t' << mapID << '\t' << info_length << '\t' << (U32)info_value_type << std::endl;

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < info_length; ++j){
				const S32 t = entry.getInteger(info_value_type, internal_pos);
				target_container += t;
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
		const U32 mapID = this->format_fields.setGet(val);

		if(this->format_containers[mapID].n_entries == 0){
			this->format_containers[mapID].setStrideSize(info_length);
		}

		if(mapID == 0){
			switch(info_value_type){
			case 1: internal_pos += this->header_->samples * sizeof(SBYTE) * info_length; break;
			case 2: internal_pos += this->header_->samples * sizeof(S16)   * info_length; break;
			case 3: internal_pos += this->header_->samples * sizeof(S32)   * info_length; break;
			}
			continue;
		}

		stream_container& target_container = this->format_containers[mapID];
		if(this->format_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 1) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 2) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 3) target_container.setType(Core::CORE_TYPE::TYPE_32B);
			else if(info_value_type == 5) target_container.setType(Core::CORE_TYPE::TYPE_FLOAT);
			else if(info_value_type == 7) target_container.setType(Core::CORE_TYPE::TYPE_8B);
			if(info_value_type != 5) target_container.header.controller.signedness = 1;
		}
		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		if(info_value_type <= 3){
			for(U32 j = 0; j < this->header_->samples*info_length; ++j){
				target_container += entry.getInteger(info_value_type, internal_pos);
			}
		}
		// Floats
		else if(info_value_type == 5){
			for(U32 j = 0; j < this->header_->samples*info_length; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(info_value_type == 7){
			for(U32 j = 0; j < this->header_->samples*info_length; ++j){
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
	const U64 hash_filter_vector = entry.hashFilter();

	U32 mapID = 0;
	if(this->format_patterns.getRaw(hash_filter_vector, mapID)){

	} else {
		std::vector<U32> ret_pattern;
		for(U32 i = 0; i < entry.filterPointer; ++i)
			ret_pattern.push_back(entry.filterID[i]);

		mapID = this->format_patterns.size();
		assert(mapID < 65536);
		this->format_patterns.set(ret_pattern, hash_filter_vector);
	}

	// Store this map in the meta
	meta.FILTER_map_ID = mapID;



	// Hash INFO pattern
	const U64 hash_info_vector = entry.hashInfo();

	mapID = 0;
	if(this->info_patterns.getRaw(hash_info_vector, mapID)){

	} else {
		std::vector<U32> ret_pattern;
		for(U32 i = 0; i < entry.infoPointer; ++i)
			ret_pattern.push_back(entry.infoID[i]);

		mapID = this->info_patterns.size();
		assert(mapID < 65536);
		this->info_patterns.set(ret_pattern, hash_info_vector);
	}

	// Store this map in the meta
	meta.INFO_map_ID = mapID;


	// Hash FORMAT pattern
	const U64 hash_format_vector = entry.hashFormat();

	 mapID = 0;
	if(this->format_patterns.getRaw(hash_format_vector, mapID)){
	} else {
		std::vector<U32> ret_pattern;
		for(U32 i = 0; i < entry.formatPointer; ++i)
			ret_pattern.push_back(entry.formatID[i]);

		mapID = this->format_patterns.size();
		assert(mapID < 65536);
		this->format_patterns.set(ret_pattern, hash_format_vector);
	}

	// Store this map in the meta
	meta.FORMAT_map_ID = mapID;

	// Return
	return true;
}

} /* namespace Tomahawk */
