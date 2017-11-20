#include <fstream>
#include "Importer.h"

namespace Tachyon {

Importer::Importer(std::string inputFile, std::string outputPrefix, const U32 checkpoint) :
	checkpoint_size(checkpoint),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader_(inputFile),
	writer_(),
	header_(nullptr)
{
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
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tachyon::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
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
	this->block.resize(resize_to);


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

		this->block.index_entry.controller.hasGT = 1;
		this->block.index_entry.controller.isDiploid = 1;

		std::cerr <<"META-HOT\t" << this->block.meta_hot_container.buffer_data.size() << '\t' << 0 << std::endl;
		std::cerr <<"META-COLD\t" << this->block.meta_cold_container.buffer_data.size() << '\t' << 0 << std::endl;
		std::cerr <<"GT RLE\t" << this->block.gt_rle_container.buffer_data.size() << '\t' << 0 << std::endl;
		std::cerr <<"GT SIMPLE\t" << this->block.gt_simple_container.buffer_data.size() << '\t' << 0 << std::endl;

		// Update head meta
		this->block.index_entry.n_info_streams = this->info_fields.size();
		this->block.index_entry.n_format_streams = this->format_fields.size();
		this->block.index_entry.n_filter_streams = this->filter_fields.size();
		this->block.index_entry.n_variants = reader.size();
		this->block.allocateOffsets(this->info_fields.size(), this->format_fields.size(), this->filter_fields.size());

		Index::IndexBlockEntryBase testBase;
		std::cerr << "external: " << testBase.contigID << std::endl;
		Index::IndexBlockEntry test;
		std::cerr << "external: " << test.contigID << std::endl;
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_INFO,   this->info_fields,   this->info_patterns);
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_FORMAT, this->format_fields, this->format_patterns);
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_FILTER, this->filter_fields, this->filter_patterns);

		std::cerr << "PATTERNS: " << this->info_patterns.size() << '\t' << this->format_patterns.size() << '\t' << this->filter_patterns.size() << std::endl;
		std::cerr << "VALUES: " << this->info_fields.size() << '\t' << this->format_fields.size() << '\t' << this->filter_fields.size() << std::endl;
		std::cerr << "INFO: " << std::endl;
		this->block.updateContainer(Index::IndexBlockEntry::INDEX_INFO,   this->info_fields);
		this->block.updateContainer(Index::IndexBlockEntry::INDEX_FORMAT, this->format_fields);
		this->block.updateContainer(Index::IndexBlockEntry::INDEX_FILTER, this->filter_fields);

		const size_t curPos = this->writer_.streamTomahawk.tellp();
		this->writer_.streamTomahawk << this->block;
		std::cerr << "Header size: " << (size_t)this->writer_.streamTomahawk.tellp() - curPos << std::endl;
		std::cerr << "cid: " << this->block.index_entry.contigID << '\t' << this->block.index_entry.info_offsets[0].key << '\t' << this->block.index_entry.info_offsets[0].header.uLength << '\t' << &this->block.index_entry.info_offsets[0] << std::endl;

		this->permutator.reset();
		this->resetHashes();
		this->block.clear();
		std::cerr << "end of reset" << std::endl;
	}

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
	if(!this->encoder.Encode(entry, meta, this->block.gt_rle_container, this->block.gt_simple_container, n_runs, this->permutator.manager.get()))
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
	case 1: this->block.meta_hot_container += (BYTE)n_runs; break;
	case 2: this->block.meta_hot_container += (U16)n_runs; break;
	case 4: this->block.meta_hot_container += (U32)n_runs; break;
	case 8: this->block.meta_hot_container += (U64)n_runs; break;
	default:
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word-size!" << std::endl;
		exit(1); // unrecoverable error
	}

	// Complex meta data
	Core::EntryColdMeta test;
	if(!test.write(entry, this->block.meta_cold_container)){
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
		this->filter_fields.setGet(val);
	}

	// At INFO
	U32 info_length;
	BYTE info_value_type;
	while(entry.nextInfo(val, info_length, info_value_type, internal_pos)){
		// Hash INFO values
		const U32 mapID = this->info_fields.setGet(val);

		//
		stream_container& target_container = this->block.info_containers[mapID];
		if(this->block.info_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			target_container.header_stride.controller.type = Core::CORE_TYPE::TYPE_32B;
			target_container.header_stride.controller.signedness = 0;
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
		// These are BCF value types
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
		const U32 mapID = this->format_fields.setGet(val);

		//
		stream_container& target_container = this->block.format_containers[mapID];
		if(this->block.format_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			target_container.header_stride.controller.type = Core::CORE_TYPE::TYPE_32B;
			target_container.header_stride.controller.signedness = 0;
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
		// These are BCF value types
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

	// Hash FILTER pattern
	const U64 hash_filter_vector = entry.hashFilter();

	U32 mapID = 0;
	if(this->filter_patterns.getRaw(hash_filter_vector, mapID)){

	} else {
		std::vector<U32> ret_pattern;
		for(U32 i = 0; i < entry.filterPointer; ++i)
			ret_pattern.push_back(entry.filterID[i]);

		mapID = this->filter_patterns.size();
		assert(mapID < 65536);
		this->filter_patterns.set(ret_pattern, hash_filter_vector);
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
