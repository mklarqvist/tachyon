#include <fstream>

#include "../algorithm/digital_digest.h"
#include "meta_cold.h"
#include "footer/footer.h"
#include "../containers/checksum_container.h"
#include "../algorithm/timer.h"
#include "variant_importer.h"

namespace tachyon {

#define IMPORT_ASSERT 1

VariantImporter::VariantImporter(std::string inputFile,
		           std::string outputPrefix,
                     const U32 checkpoint_n_snps,
                  const double checkpoint_bases) :
	GT_available_(false),
	permute(true),
	checkpoint_n_snps(checkpoint_n_snps),
	checkpoint_bases(checkpoint_bases),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	writer(),
	header(nullptr)
{
}

VariantImporter::~VariantImporter(){ }

void VariantImporter::resetHashes(void){
	this->info_fields.clear();
	this->info_patterns.clear();
	this->format_fields.clear();
	this->format_patterns.clear();
	this->filter_fields.clear();
	this->filter_patterns.clear();
}

bool VariantImporter::Build(){
	std::ifstream temp(this->inputFile, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << utility::timestamp("ERROR", "IMPORT")  << "Failed to open file (" << this->inputFile << ")..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if((BYTE)tempData[0] == io::constants::GZIP_ID1 && (BYTE)tempData[1] == io::constants::GZIP_ID2){
		if(!this->BuildBCF()){
			std::cerr << utility::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << utility::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool VariantImporter::BuildBCF(void){
	bcf_reader_type reader;
	if(!reader.open(this->inputFile)){
		std::cerr << utility::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	this->header = &reader.header;

	// Spawn RLE controller and update PPA controller
	this->encoder.setSamples(this->header->samples);
	this->block.ppa_manager.setSamples(this->header->samples);
	this->permutator.manager = &this->block.ppa_manager;
	this->permutator.setSamples(this->header->samples);

	if(!this->writer.open(this->outputPrefix)){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	core::VariantHeader header(*this->header);
	header.write(this->writer.stream);
	this->GT_available_ = header.has_format_field("GT");
	for(U32 i = 0; i < this->header->format_map.size(); ++i){
		if(this->header->format_map[i].ID == "GT"){
			reader.map_gt_id = this->header->format_map[i].IDX;
		}
	}

	this->block.header.controller.hasGT = this->GT_available_;

	// Resize containers
	const U32 resize_to = this->checkpoint_n_snps * sizeof(U32) * this->header->samples * 100;
	this->block.resize(resize_to);
	if(this->header->samples == 0){
		this->block.resize(this->checkpoint_n_snps * sizeof(U32) * 100);
	}

	// Digest controller
	containers::ChecksumContainer checksums;
	if(checksums.allocate(this->header->info_map.size() + this->header->format_map.size() + this->header->filter_map.size()) == false){
		std::cerr << "failed to allocate" << std::endl;
		return false;
	}
	//tachyon::DigitalDigestPair* digests = new tachyon::DigitalDigestPair[this->header->map.size()];

	// Start import
	U32 previousFirst    = 0;
	U32 previousLast     = 0;
	S32 previousContigID = -1;
	U64 n_variants_read  = 0;

	// Index
	index_entry_type current_index_entry;

	// Begin import
	// Get BCF entries
	algorithm::Timer timer; timer.Start();
	while(true){
		if(!reader.getVariants(this->checkpoint_n_snps, this->checkpoint_bases)){
			break;
		}

		// Debug assertion
#if IMPORT_ASSERT == 1
		if(reader.front().body->CHROM == previousContigID){
			if(!(reader.front().body->POS >= previousFirst && reader.front().body->POS >= previousLast)){
				std::cerr << utility::timestamp("ERROR","IMPORT") << reader.front().body->POS << '/' << previousFirst << '/' << previousLast << std::endl;
				std::cerr << reader[reader.n_entries].body->POS << std::endl;
				exit(1);
			}
		}
#endif
		this->block.header.contigID    = reader.front().body->CHROM;
		this->block.header.minPosition = reader.front().body->POS;
		this->block.header.maxPosition = reader.back().body->POS;
		this->block.header.controller.hasGT = this->GT_available_;
		this->block.header.controller.hasGTPermuted = this->permute;
		// if there is 0 or 1 samples then GT data is never permuted
		if(header.getSampleNumber() <= 1)
			this->block.header.controller.hasGTPermuted = false;

		// Permute GT if GT is available and the appropriate flag is triggered
		if(this->block.header.controller.hasGT && this->block.header.controller.hasGTPermuted){
			if(!this->permutator.build(reader)){
				std::cerr << utility::timestamp("ERROR","PERMUTE") << "Failed to complete..." << std::endl;
				return false;
			}
		}

		// Perform parsing of BCF entries in memory
		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->add(reader[i])){
				std::cerr << utility::timestamp("ERROR","IMPORT") << "Failed to add BCF entry..." << std::endl;
				return false;
			}
		}

		// Stats
		n_variants_read += reader.size();

		// Update head meta
		this->block.header.controller.hasGT = this->GT_available_;
		this->block.footer.n_info_streams   = this->info_fields.size();
		this->block.footer.n_filter_streams = this->filter_fields.size();
		this->block.footer.n_format_streams = this->format_fields.size();
		this->block.header.n_variants       = reader.size();
		this->block.allocateDiskOffsets(this->info_fields.size(), this->format_fields.size(), this->filter_fields.size());
		this->block.updateBaseContainers();
		this->block.updateContainerSet(containers::DataBlockFooter::INDEX_INFO);
		this->block.updateContainerSet(containers::DataBlockFooter::INDEX_FORMAT);

		// Perform compression using standard parameters
		if(!this->compression_manager.compress(this->block)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Digests
		//if(checksums.update(this->block, this->header->mapTable) == false){
		//	std::cerr << utility::timestamp("ERROR","CHECKSUM") << "Failed to update!" << std::endl;
		//	return false;
		//}

		// Todo: abstraction
		// Perform writing and update index
		this->block.footer.constructBitVector(containers::DataBlockFooter::INDEX_INFO,   this->info_fields,   this->info_patterns);
		this->block.footer.constructBitVector(containers::DataBlockFooter::INDEX_FILTER, this->filter_fields, this->filter_patterns);
		this->block.footer.constructBitVector(containers::DataBlockFooter::INDEX_FORMAT, this->format_fields, this->format_patterns);

		current_index_entry.byte_offset     = this->writer.stream.tellp();

		this->block.write(this->writer.stream, this->stats_basic, this->stats_info, this->stats_format);
		current_index_entry.byte_offset_end = this->writer.stream.tellp();
		current_index_entry.contigID        = reader.front().body->CHROM;
		current_index_entry.minPosition     = reader.front().body->POS;
		current_index_entry.maxPosition     = reader.back().body->POS;
		current_index_entry.n_variants      = reader.size();
		this->writer.index += current_index_entry;
		current_index_entry.reset();
		++this->writer.n_blocks_written;
		this->writer.n_variants_written += reader.size();

		if(!SILENT){
			std::cerr << utility::timestamp("PROGRESS") <<
			std::setfill(' ') << std::setw(10) << this->writer.n_variants_written << ' ' <<
			std::setfill(' ') << std::setw(5) << reader.size() << '\t' <<
			header.contigs[reader.front().body->CHROM].name << ":" << reader.front().body->POS+1 << "->" << reader.back().body->POS+1 << "\t" <<
			std::setfill(' ') << std::setw(8) << (double)reader.stream.tellg()/reader.filesize*100 << "%" << ' ' <<
			timer.ElapsedString() << std::endl;
		}

		// Reset and update
		this->resetHashes();
		this->block.clear();
		this->permutator.reset();
		this->writer.stream.flush();
		previousContigID = reader.front().body->CHROM;
		previousFirst    = reader.front().body->POS;
		previousLast     = reader.back().body->POS;
	}
	// Done importing
	this->writer.stream.flush();

	core::Footer footer;
	footer.offset_end_of_data = this->writer.stream.tellp();
	footer.n_blocks           = this->writer.n_blocks_written;
	footer.n_variants         = this->writer.n_variants_written;

	// Write index
	this->writer.WriteIndex();
	const U64 index_ends = this->writer.stream.tellp();

	// Finalize SHA-512 digests
	// Write digests
	checksums.finalize();
	this->writer.stream << checksums;

	this->writer.stream.flush();
	const U64 digests_ends = this->writer.stream.tellp();

	this->writer.stream << footer;
	this->writer.stream.flush();

	if(!SILENT){
		std::cerr << "Header:    " << utility::ToPrettyString(this->stats_basic[0].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[0].cost_uncompressed) << '\t' << (double)this->stats_basic[0].cost_uncompressed/this->stats_basic[0].cost_compressed << "-fold" << std::endl;
		std::cerr << "PPA:       " << utility::ToPrettyString(this->stats_basic[1].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[1].cost_uncompressed) << '\t' << (double)this->stats_basic[1].cost_uncompressed/this->stats_basic[1].cost_compressed << "-fold" << std::endl;
		std::cerr << "Meta hot:  " << utility::ToPrettyString(this->stats_basic[2].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[2].cost_uncompressed) << '\t' << (double)this->stats_basic[2].cost_uncompressed/this->stats_basic[2].cost_compressed << "-fold" << std::endl;
		std::cerr << "Meta cold: " << utility::ToPrettyString(this->stats_basic[3].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[3].cost_uncompressed) << '\t' << (double)this->stats_basic[3].cost_uncompressed/this->stats_basic[3].cost_compressed << "-fold" << std::endl;
		std::cerr << "GT RLE:    " << utility::ToPrettyString(this->stats_basic[4].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[4].cost_uncompressed) << '\t' << (double)this->stats_basic[4].cost_uncompressed/this->stats_basic[4].cost_compressed << "-fold" << std::endl;
		std::cerr << "GT Packed: " << utility::ToPrettyString(this->stats_basic[5].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[5].cost_uncompressed) << '\t' << (double)this->stats_basic[5].cost_uncompressed/this->stats_basic[5].cost_compressed << "-fold" << std::endl;
		std::cerr << "GT Support:" << utility::ToPrettyString(this->stats_basic[6].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[6].cost_uncompressed) << '\t' << (double)this->stats_basic[6].cost_uncompressed/this->stats_basic[6].cost_compressed << "-fold" << std::endl;
		std::cerr << "Meta IDs:  " << utility::ToPrettyString(this->stats_basic[7].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[7].cost_uncompressed) << '\t' << (double)this->stats_basic[7].cost_uncompressed/this->stats_basic[7].cost_compressed << "-fold" << std::endl;
		std::cerr << "INFO:      " << utility::ToPrettyString(this->stats_basic[8].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[8].cost_uncompressed) << '\t' << (double)this->stats_basic[8].cost_uncompressed/this->stats_basic[8].cost_compressed << "-fold" << std::endl;
		std::cerr << "FORMAT:    " << utility::ToPrettyString(this->stats_basic[9].cost_compressed) << "\t" << utility::ToPrettyString(this->stats_basic[9].cost_uncompressed) << '\t' << (double)this->stats_basic[9].cost_uncompressed/this->stats_basic[9].cost_compressed << "-fold" << std::endl;
		std::cerr << "Checksums: " << utility::ToPrettyString(digests_ends - index_ends) << std::endl;

		for(U32 i = 0; i < header.header_magic.n_info_values; ++i){
			const U64 cost = this->stats_info[i].cost_uncompressed;
			if(cost)
				std::cerr << std::setw(14) << header.info_fields[i].ID << "\t" << std::setw(14) << utility::ToPrettyString(this->stats_info[i].cost_compressed) << "\t" << std::setw(14) << utility::ToPrettyString(this->stats_info[i].cost_uncompressed) << '\t' << std::setw(14) << (double)this->stats_info[i].cost_uncompressed/this->stats_info[i].cost_compressed << "-fold" << std::endl;
			else
				std::cerr << std::setw(14) << header.info_fields[i].ID << '\t' << "unused" << std::endl;
		}

		for(U32 i = 0; i < header.header_magic.n_format_values; ++i){
			const U64 cost = this->stats_format[i].cost_uncompressed;
			if(cost)
				std::cerr << std::setw(14) << header.format_fields[i].ID << '\t' << std::setw(14) << utility::ToPrettyString(this->stats_format[i].cost_compressed) << "\t" << std::setw(14) << utility::ToPrettyString(this->stats_format[i].cost_uncompressed) << '\t' << std::setw(14) << (double)this->stats_format[i].cost_uncompressed/this->stats_format[i].cost_compressed << "-fold" << std::endl;
			else
				std::cerr << std::setw(14) << header.format_fields[i].ID << '\t' << "unused" << std::endl;
		}

		std::cerr << utility::timestamp("PROGRESS") << "Wrote: " << utility::ToPrettyString(this->writer.n_variants_written) << " variants in " << utility::ToPrettyString(this->writer.n_blocks_written) << " blocks in " << timer.ElapsedString() << " to " << utility::toPrettyDiskString((U64)this->writer.stream.tellp()) << std::endl;
	}

	// All done
	return(true);
}

bool VariantImporter::add(bcf_entry_type& entry){
	// Assert position is in range
	if(entry.body->POS + 1 > this->header->getContig(entry.body->CHROM).bp_length){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << this->header->getContig(entry.body->CHROM).name << ':' << entry.body->POS+1 << " > reported max size of contig (" << this->header->getContig(entry.body->CHROM).bp_length << ")..." << std::endl;
		return false;
	}

	meta_type meta;
	meta.position = entry.body->POS;
	meta.contigID = entry.body->CHROM;
	meta.ref_alt  = entry.ref_alt;
	meta.controller.simple_snv = entry.isSimple();

	// GT encoding if available
	if(entry.hasGenotypes){
		meta.controller.gt_available = true;
		if(!this->encoder.Encode(entry,
								 meta,
								 this->block.gt_rle_container,
								 this->block.gt_simple_container,
								 this->block.gt_support_data_container,
								 this->permutator.manager->get()))
		{
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode GT..." << std::endl;
			return false;
		}
	} else {
		meta.controller.gt_available = false;
	}

	if(!this->parseBCFBody(meta, entry)){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode BCF body..." << std::endl;
		return false;
	}

	// Complex meta data
	core::MetaCold test;
	if(!test.write(entry, this->block.meta_cold_container)){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to write complex meta!" << std::endl;
		return false;
	}

	++this->block.meta_cold_container.n_entries;
	++this->block.meta_cold_container.n_additions;
	++this->block.meta_hot_container.n_entries;
	++this->block.meta_hot_container.n_additions;

	// Update number of entries in block
	++this->index_entry.n_variants;

	// Push meta
	this->block.meta_hot_container.buffer_data_uncompressed += meta;

	return true;
}

bool VariantImporter::parseBCFBody(meta_type& meta, bcf_entry_type& entry){
	for(U32 i = 0; i < entry.filterPointer; ++i){
		assert(entry.filterID[i].mapID != -1);
		this->filter_fields.setGet(this->header->filter_remap[entry.filterID[i].mapID]);
	}

	for(U32 i = 0; i < entry.infoPointer; ++i){
		assert(entry.infoID[i].mapID != -1);
		const U32 mapID = this->info_fields.setGet(this->header->info_remap[entry.infoID[i].mapID]);

		stream_container& target_container = this->block.info_containers[mapID];
		if(this->block.info_containers[mapID].size() == 0){
			target_container.setStrideSize(entry.infoID[i].l_stride);
			target_container.header.stride_header.controller.type       = YON_TYPE_32B;
			target_container.header.stride_header.controller.signedness = 0;
			// Set all integer types to U32
			// Change to smaller type later if required
			if(entry.infoID[i].primitive_type == 0)      target_container.setType(YON_TYPE_32B);
			else if(entry.infoID[i].primitive_type == 1) target_container.setType(YON_TYPE_32B);
			else if(entry.infoID[i].primitive_type == 2) target_container.setType(YON_TYPE_32B);
			else if(entry.infoID[i].primitive_type == 3) target_container.setType(YON_TYPE_32B);
			else if(entry.infoID[i].primitive_type == 5) target_container.setType(YON_TYPE_FLOAT);
			else if(entry.infoID[i].primitive_type == 7) target_container.setType(YON_TYPE_CHAR);
			else {
				std::cerr << "not possible" << std::endl;
				exit(1);
			}
			if(entry.infoID[i].primitive_type != 5)      target_container.header.data_header.controller.signedness = 1;
		}

		++target_container;
		if(!target_container.checkStrideSize(entry.infoID[i].l_stride))
			target_container.triggerMixedStride();

		target_container.addStride(entry.infoID[i].l_stride);

		// Flags and integers
		// These are BCF value types
		U32 internal_pos = entry.infoID[i].l_offset;
		if(entry.infoID[i].primitive_type <= 3){
			for(U32 j = 0; j < entry.infoID[i].l_stride; ++j){
				target_container += entry.getInteger(entry.infoID[i].primitive_type, internal_pos);
			}
		}
		// Floats
		else if(entry.infoID[i].primitive_type == 5){
			for(U32 j = 0; j < entry.infoID[i].l_stride; ++j){
				target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(entry.infoID[i].primitive_type == 7){
			for(U32 j = 0; j < entry.infoID[i].l_stride; ++j){
				target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible in info: " << (int)entry.infoID[i].primitive_type << std::endl;
			exit(1);
		}
	}

	for(U32 i = 0; i < entry.formatPointer; ++i){
		assert(entry.formatID[i].mapID != -1);

		const U32 mapID = this->format_fields.setGet(this->header->format_remap[entry.formatID[i].mapID]);
		U32 internal_pos = entry.formatID[i].l_offset;

		// First value is always genotypes if there are any
		if(entry.hasGenotypes == true && i == 0)
			continue;

		// Hash INFO values
		stream_container& target_container = this->block.format_containers[mapID];
		if(this->block.format_containers[mapID].size() == 0){
			target_container.setStrideSize(entry.formatID[i].l_stride);
			target_container.header.stride_header.controller.type       = YON_TYPE_32B;
			target_container.header.stride_header.controller.signedness = 0;
			// Set all integer types to U32
			// Change to smaller type later if required
			if(entry.formatID[i].primitive_type == 0)      target_container.setType(YON_TYPE_32B);
			else if(entry.formatID[i].primitive_type == 1) target_container.setType(YON_TYPE_32B);
			else if(entry.formatID[i].primitive_type == 2) target_container.setType(YON_TYPE_32B);
			else if(entry.formatID[i].primitive_type == 3) target_container.setType(YON_TYPE_32B);
			else if(entry.formatID[i].primitive_type == 5) target_container.setType(YON_TYPE_FLOAT);
			else if(entry.formatID[i].primitive_type == 7) target_container.setType(YON_TYPE_CHAR);
			else {
				std::cerr << "not possible" << std::endl;
				exit(1);
			}
			if(entry.formatID[i].primitive_type != 5)      target_container.header.data_header.controller.signedness = 1;
		}

		++target_container;
		if(!target_container.checkStrideSize(entry.formatID[i].l_stride))
			target_container.triggerMixedStride();

		target_container.addStride(entry.formatID[i].l_stride);

		// Flags and integers
		// These are BCF value types
		if(entry.formatID[i].primitive_type <= 3){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container += entry.getInteger(entry.formatID[i].primitive_type, internal_pos);
			}
		}
		// Floats
		else if(entry.formatID[i].primitive_type == 5){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container += entry.getFloat(internal_pos);
			}
		}
		// Chars
		else if(entry.formatID[i].primitive_type == 7){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container += entry.getChar(internal_pos);
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)entry.formatID[i].primitive_type << std::endl;
			std::cerr << utility::timestamp("LOG") << entry.formatID[i].mapID << '\t' << entry.formatID[i].l_stride << '\t' << (int)entry.formatID[i].primitive_type << '\t' << internal_pos << '/' << entry.l_data << std::endl;
			exit(1);
		}
	}

	if(entry.filterPointer){
		// Hash FILTER pattern
		const U64 hash_filter_vector = entry.hashFilter();

		U32 mapID = 0;
		if(this->filter_patterns.getRaw(hash_filter_vector, mapID)){

		} else {
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.filterPointer; ++i)
				ret_pattern.push_back(this->header->filter_remap[entry.filterID[i].mapID]);

			mapID = this->filter_patterns.size();
			assert(mapID < 65536);
			if(!this->filter_patterns.set(ret_pattern, hash_filter_vector)){
				std::cerr << "failed to insert filter: " << ret_pattern.size() << " and " << hash_filter_vector << std::endl;
				std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


				for(size_t i = 0; i < ret_pattern.size(); ++i){
					std::cerr << ret_pattern[i] << std::endl;
				}
				exit(1);
			}
		}

		// Store this map in the meta
		this->block.meta_filter_map_ids += (S32)mapID;
		++this->block.meta_filter_map_ids;
	} else {
		this->block.meta_filter_map_ids += -1;
		++this->block.meta_filter_map_ids;
	}

	if(entry.infoPointer){
		U32 mapID = 0;
		// Hash INFO pattern

		const U64 hash_info_vector = entry.hashInfo();
		//const U64 hash_info_vector = 0;

		if(this->info_patterns.getRaw(hash_info_vector, mapID)){

		} else {
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.infoPointer; ++i)
				ret_pattern.push_back(this->header->info_remap[entry.infoID[i].mapID]);

			mapID = this->info_patterns.size();
			assert(mapID < 65536);
			if(!this->info_patterns.set(ret_pattern, hash_info_vector)){
				std::cerr << "failed to insert info: " << ret_pattern.size() << " and " << hash_info_vector << std::endl;
				std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


				for(size_t i = 0; i < ret_pattern.size(); ++i){
					std::cerr << ret_pattern[i] << std::endl;
				}
				exit(1);
			}
		}

		// Store this map in the meta
		//meta.INFO_map_ID = mapID;
		this->block.meta_info_map_ids += (S32)mapID;
		++this->block.meta_info_map_ids;
	} else {
		this->block.meta_info_map_ids += -1;
		++this->block.meta_info_map_ids;
	}

	if(entry.formatPointer){
		U32 mapID = 0;
		// Hash FORMAT pattern
		const U64 hash_format_vector = entry.hashFormat();

		if(this->format_patterns.getRaw(hash_format_vector, mapID)){

		} else {
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.formatPointer; ++i)
				ret_pattern.push_back(this->header->format_remap[entry.formatID[i].mapID]);

			mapID = this->format_patterns.size();
			assert(mapID < 65536);
			if(!this->format_patterns.set(ret_pattern, hash_format_vector)){
				std::cerr << "failed to insert format: " << ret_pattern.size() << " and " << hash_format_vector << std::endl;
				std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


				for(size_t i = 0; i < ret_pattern.size(); ++i){
					std::cerr << ret_pattern[i] << std::endl;
				}
				exit(1);
			}
		}

		// Store this map in the meta
		//meta.FORMAT_map_ID = mapID;
		this->block.meta_format_map_ids += (S32)mapID;
		++this->block.meta_format_map_ids;
	} else {
		this->block.meta_format_map_ids += -1;
		++this->block.meta_format_map_ids;
	}

	// Return
	return true;
}

} /* namespace Tachyon */
