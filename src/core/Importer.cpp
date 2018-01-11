#include <fstream>

#include "Importer.h"
#include "base/MetaCold.h"
#include "../algorithm/DigitalDigestController.h"

namespace Tachyon {

#define IMPORT_ASSERT 1

Importer::Importer(std::string inputFile,
		           std::string outputPrefix,
                     const U32 checkpoint_n_snps,
                  const double checkpoint_bases) :
	permute(true),
	checkpoint_n_snps(checkpoint_n_snps),
	checkpoint_bases(checkpoint_bases),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader(inputFile),
	writer(),
	header(nullptr)
{
}

Importer::~Importer(){ }

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

	this->header = &reader.header;
	if(this->header->samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "No samples detected in header..." << std::endl;
		return false;
	}

	if(this->header->samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tachyon::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller and update PPA controller
	this->encoder.setSamples(this->header->samples);
	this->block.ppa_manager.setSamples(this->header->samples);
	this->permutator.manager = &this->block.ppa_manager;
	this->permutator.setSamples(this->header->samples);

	if(!this->writer.Open(this->outputPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	if(!this->writer.WriteHeader()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to write header..." << std::endl;
		return false;
	}
	this->writer.stream << *this->header;

	// Resize containers
	const U32 resize_to = this->checkpoint_n_snps * sizeof(U32) * this->header->samples * 100;
	this->block.resize(resize_to);

	// Digest controller
	Algorithm::DigitalDigestController* digests = new Algorithm::DigitalDigestController[this->header->map.size()];
	for(U32 i = 0; i < this->header->map.size(); ++i){
		if(!digests[i].initialize()){
			std::cerr << "failed to init sha512" << std::endl;
			return false;
		}
	}

	// Start import
	U32 previousFirst    = 0;
	U32 previousLast     = 0;
	S32 previousContigID = -1;
	U64 n_variants_read  = 0;

	// Index
	Index::IndexEntry current_index_entry;

	// Begin import
	// Get BCF entries
	bcf_entry_type t;
	while(true){
		if(!reader.getVariants(this->checkpoint_n_snps, this->checkpoint_bases)){
			break;
		}

		// Debug assertion
#if IMPORT_ASSERT == 1
		if(reader.first().body->CHROM == previousContigID){
			if(!(reader.first().body->POS >= previousFirst && reader.first().body->POS >= previousLast)){
				std::cerr << Helpers::timestamp("ERROR","IMPORT") << reader.first().body->POS << '/' << previousFirst << '/' << previousLast << std::endl;
				std::cerr << reader[reader.n_entries].body->POS << std::endl;
				exit(1);
			}
		}
#endif
		std::cerr << Helpers::timestamp("DEBUG") << "n_variants: " << reader.size() << '\t' << reader.first().body->POS+1 << "->" << reader.last().body->POS+1 << std::endl;
		this->block.index_entry.contigID    = reader.first().body->CHROM;
		this->block.index_entry.minPosition = reader.first().body->POS;
		this->block.index_entry.maxPosition = reader.last().body->POS;
		this->block.index_entry.controller.hasGTPermuted = this->permute;

		// Permute or not?
		if(this->block.index_entry.controller.hasGTPermuted){
			if(!this->permutator.build(reader)){
				std::cerr << Helpers::timestamp("ERROR","PERMUTE") << "Failed to complete..." << std::endl;
				return false;
			}
		}

		// Perform parsing of BCF entries in memory
		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->parseBCFLine(reader[i])){
				std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Failed to parse BCF body..." << std::endl;
				return false;
			}
		}

		// Stats
		n_variants_read += reader.size();

		assert(this->block.gt_support_data_container.n_entries == reader.size());
		assert(this->block.meta_info_map_ids.n_entries         == reader.size());
		assert(this->block.meta_filter_map_ids.n_entries       == reader.size());
		assert(this->block.meta_format_map_ids.n_entries       == reader.size());

		// Update head meta
		this->block.index_entry.controller.hasGT     = true;
		this->block.index_entry.controller.isDiploid = true;
		this->block.index_entry.n_info_streams       = this->info_fields.size();
		this->block.index_entry.n_filter_streams     = this->filter_fields.size();
		this->block.index_entry.n_format_streams     = this->format_fields.size();
		this->block.index_entry.n_variants           = reader.size();
		this->block.allocateOffsets(this->info_fields.size(), this->format_fields.size(), this->filter_fields.size());
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_INFO,   this->info_fields,   this->info_patterns);
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_FILTER, this->filter_fields, this->filter_patterns);
		this->block.index_entry.constructBitVector(Index::IndexBlockEntry::INDEX_FORMAT, this->format_fields, this->format_patterns);
		this->block.updateBaseContainers();
		this->block.updateContainerSet(Index::IndexBlockEntry::INDEX_INFO);
		this->block.updateContainerSet(Index::IndexBlockEntry::INDEX_FORMAT);

		// Todo: abstraction
		// Digests
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			if(!digests[this->header->mapTable[this->block.index_entry.info_offsets[i].key]].updateUncompressed(this->block.info_containers[i])){
				std::cerr << Helpers::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
				return false;
			}
		}

		for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
			if(!digests[this->header->mapTable[this->block.index_entry.format_offsets[i].key]].updateUncompressed(this->block.format_containers[i])){
				std::cerr << Helpers::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
				return false;
			}
		}

		// Perform compression using standard parameters
		if(!compression_manager.compress(this->block)){
			std::cerr << Helpers::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Todo: abstraction
		// Perform writing and update index
		this->block.updateOffsets();
		current_index_entry.byte_offset = this->writer.stream.tellp();
		this->block.write(this->writer.stream, this->import_compressed_stats, this->import_uncompressed_stats);
		current_index_entry.byte_offset_end = this->writer.stream.tellp();
		current_index_entry.contigID    = reader.first().body->CHROM;
		current_index_entry.minPosition = reader.first().body->POS;
		current_index_entry.maxPosition = reader.last().body->POS;
		current_index_entry.n_variants  = reader.size();
		this->writer.index += current_index_entry;
		current_index_entry.reset();

		// Reset and update
		this->resetHashes();
		this->block.clear();
		this->permutator.reset();
		this->writer.stream.flush();
		previousContigID = reader.first().body->CHROM;
		previousFirst    = reader.first().body->POS;
		previousLast     = reader.last().body->POS;
	}
	// Done importing
	this->writer.stream.flush();
	const U64 data_ends = this->writer.stream.tellp();

	// Write index
	this->writer.WriteIndex();
	const U64 index_ends = this->writer.stream.tellp();

	// Finalize SHA-512 digests
	const U64 digests_start = this->writer.stream.tellp();

	for(U32 i = 0; i < this->header->map.size(); ++i){
		digests[i].finalize();
		this->writer.stream << digests[i];
		//std::cerr << std::hex;
		//for(U32 j = 0; j < 64; ++j)
		//	std::cerr << std::hex << (int)digests[i].sha512_digest[j];

		//std::cerr << std::dec << std::endl;
	}
	delete [] digests;
	this->writer.stream.flush();
	const U64 digests_ends = this->writer.stream.tellp();


	// Place markers
	this->writer.stream.write(reinterpret_cast<const char* const>(&digests_start), sizeof(U64));
	this->writer.WriteFinal(data_ends);

	std::cout
	    << "Header:    " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_header_cost) << '\n'
		<< "GT:        " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_gt_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_gt_cost) << '\n'
		<< "PPA:       " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_ppa_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_ppa_cost) << '\n'
		<< "Meta:      " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_meta_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_meta_cost) << '\n'
		<< "INFO:      " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_info_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_info_cost) << '\n'
		<< "FORMAT:    " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_format_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_format_cost) << '\n'
		<< "IDs:       " << Helpers::toPrettyDiskString(this->import_compressed_stats.total_special_cost) << '\t' << Helpers::toPrettyDiskString(this->import_uncompressed_stats.total_special_cost) << '\n'
		<< "Index:     " << Helpers::toPrettyDiskString(index_ends - data_ends) << '\n'
		<< "Checksums: " << Helpers::toPrettyDiskString(digests_ends - index_ends) << '\n'
		<< "Total:     " << Helpers::toPrettyDiskString((U64)this->writer.stream.tellp()) << std::endl;

	// All done
	return(true);
}

bool Importer::parseBCFLine(bcf_entry_type& entry){
	// Assert position is in range
	if(entry.body->POS + 1 > this->header->getContig(entry.body->CHROM).bp_length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << this->header->getContig(entry.body->CHROM).name << ':' << entry.body->POS+1 << " > reported max size of contig (" << this->header->getContig(entry.body->CHROM).bp_length << ")..." << std::endl;
		return false;
	}

	// Perform run-length encoding
	U64 n_runs = 0;

	meta_type meta;
	meta.position          = entry.body->POS - this->block.index_entry.minPosition;
	meta.ref_alt           = entry.ref_alt;

	// GT encoding
	if(!this->encoder.Encode(entry, meta, this->block.gt_rle_container, this->block.gt_simple_container, this->block.gt_support_data_container, n_runs, this->permutator.manager->get())){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to encode GT..." << std::endl;
		return false;
	}

	if(!this->parseBCFBody(meta, entry)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to encode BCF body..." << std::endl;
		return false;
	}

	// Complex meta data
	Core::MetaCold test;
	if(!test.write(entry, this->block.meta_cold_container)){
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	++this->block.meta_cold_container.n_entries;
	++this->block.meta_cold_container.n_additions;
	++this->block.meta_hot_container.n_entries;
	++this->block.meta_hot_container.n_additions;

	// Update number of entries in block
	++this->index_entry.n_variants;

	//meta.n_objects = n_runs;
	this->block.gt_support_data_container += (U32)n_runs;
	++this->block.gt_support_data_container;
	this->block.meta_hot_container.buffer_data_uncompressed += meta;

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
		//std::cerr << "FILTER: " << val << "->" << (*this->header_)[val].ID << std::endl;
	}

	// At INFO
	U32 info_length;
	BYTE info_value_type;
	while(entry.nextInfo(val, info_length, info_value_type, internal_pos)){
		// Hash INFO values
		const U32 mapID = this->info_fields.setGet(val);
		//std::cerr << Helpers::timestamp("DEBUG") << "Field: " << val << "->" << this->header_->mapTable[val] << "->" << this->header_->map[this->header_->mapTable[val]].ID << '\t' << mapID << '\t' << "length: " << info_length << '\t' <<internal_pos << "/" << (entry.body->l_shared + sizeof(U32)*2) << std::endl;

		stream_container& target_container = this->block.info_containers[mapID];
		if(this->block.info_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			target_container.header_stride.controller.type = Core::YON_TYPE_32B;
			target_container.header_stride.controller.signedness = 0;
			// Set all integer types to U32
			// Change to smaller type later if required
			if(info_value_type == 0)      target_container.setType(Core::YON_TYPE_32B);
			else if(info_value_type == 1) target_container.setType(Core::YON_TYPE_32B);
			else if(info_value_type == 2) target_container.setType(Core::YON_TYPE_32B);
			else if(info_value_type == 3) target_container.setType(Core::YON_TYPE_32B);
			else if(info_value_type == 5) target_container.setType(Core::YON_TYPE_FLOAT);
			else if(info_value_type == 7) target_container.setType(Core::YON_TYPE_CHAR);
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
			std::cerr << "impossible in info: " << (int)info_value_type << std::endl;
			exit(1);
		}
	}

#if IMPORT_ASSERT == 1
	// Assert all FILTER and INFO data have been successfully
	// parsed. This is true when the byte pointer equals the
	// start position of the FORMAT fields which are encoded
	// in the meta header structure
	assert(internal_pos == (entry.body->l_shared + sizeof(U32)*2));
#endif

	BYTE format_value_type = 0;
	while(entry.nextFormat(val, info_length, format_value_type, internal_pos)){
		const U32 mapID = this->format_fields.setGet(val);

		// If this is a GT field
		// then continue
		if(this->header->map[this->header->mapTable[val]].ID == "GT"){
			BYTE multiplier = sizeof(U32);
			switch(format_value_type){
				case(7):
				case(1): multiplier = sizeof(SBYTE); break;
				case(2): multiplier = sizeof(S16);   break;
				case(5):
				case(3): multiplier = sizeof(S32);   break;
				case(0): break; // FLAG
				default: std::cerr << "illegal" << std::endl; exit(1); return false;
			}
			internal_pos += this->header->samples * info_length * multiplier;
			continue;
		}

		// Hash INFO values
		stream_container& target_container = this->block.format_containers[mapID];
		if(this->block.format_containers[mapID].n_entries == 0){
			target_container.setStrideSize(info_length);
			target_container.header_stride.controller.type = Core::YON_TYPE_32B;
			target_container.header_stride.controller.signedness = 0;
			// Set all integer types to U32
			// Change to smaller type later if required
			if(format_value_type == 0)      target_container.setType(Core::YON_TYPE_32B);
			else if(format_value_type == 1) target_container.setType(Core::YON_TYPE_32B);
			else if(format_value_type == 2) target_container.setType(Core::YON_TYPE_32B);
			else if(format_value_type == 3) target_container.setType(Core::YON_TYPE_32B);
			else if(format_value_type == 5) target_container.setType(Core::YON_TYPE_FLOAT);
			else if(format_value_type == 7) target_container.setType(Core::YON_TYPE_CHAR);
			else {
				std::cerr << "not possible" << std::endl;
				exit(1);
			}
			if(format_value_type != 5) target_container.header.controller.signedness = 1;
		}

		++target_container;
		if(!target_container.checkStrideSize(info_length))
			target_container.setMixedStrides();

		target_container.addStride(info_length);

		// Flags and integers
		// These are BCF value types
		if(format_value_type <= 3){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					target_container += entry.getInteger(format_value_type, internal_pos);
			}
		}
		// Floats
		else if(format_value_type == 5){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					target_container += entry.getFloat(internal_pos);

			}
		}
		// Chars
		else if(format_value_type == 7){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < info_length; ++j)
					target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)format_value_type << std::endl;
			std::cerr << Helpers::timestamp("LOG") << val << '\t' << info_length << '\t' << (int)format_value_type << '\t' << internal_pos << '/' << entry.pointer << std::endl;

			exit(1);
		}
	}
#if IMPORT_ASSERT == 1
	assert(internal_pos == entry.pointer);
#endif

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
	this->block.meta_filter_map_ids += mapID;
	//meta.FILTER_map_ID = mapID;

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
	//meta.INFO_map_ID = mapID;
	this->block.meta_info_map_ids += mapID;


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
	//meta.FORMAT_map_ID = mapID;
	this->block.meta_format_map_ids += mapID;

	// Update
	++this->block.meta_info_map_ids;
	++this->block.meta_format_map_ids;
	++this->block.meta_filter_map_ids;
	this->block.meta_info_map_ids.addStride((S32)1);
	this->block.meta_format_map_ids.addStride((S32)1);
	this->block.meta_filter_map_ids.addStride((S32)1);

	// Return
	return true;
}

} /* namespace Tachyon */
