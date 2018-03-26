#include <fstream>

#include "../algorithm/digital_digest.h"
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
	if(!SILENT){
		std::cerr << utility::timestamp("PROGRESS") <<
		std::setfill(' ') << std::setw(10) << "Variants" << ' ' <<
		std::setfill(' ') << std::setw(10) << "Written" << '\t' <<
		std::setfill(' ') << std::setw(8) << "Completion" << ' ' <<
		"Elapsed " << "Contig:from->to" << std::endl;
	}

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
		this->block.header.controller.hasGT         = this->GT_available_;
		this->block.header.controller.hasGTPermuted = this->permute;
		// if there is 0 or 1 samples then GT data is never permuted
		if(header.getSampleNumber() <= 1)
			this->block.header.controller.hasGTPermuted = false;

		// test
		/*
		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->add(reader[i])){
				std::cerr << utility::timestamp("ERROR","IMPORT") << "Failed to add BCF entry..." << std::endl;
				return false;
			}

			this->permutator.update(reader[i]);
		}
		*/

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
		this->block.header.n_variants       = reader.size();
		this->block.finalize();

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

		current_index_entry.byte_offset     = this->writer.stream.tellp();
		//this->block.footer_support.buffer_data_uncompressed += this->block.footer;
		//this->compression_manager.zstd_codec.compress(this->block.footer_support);

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
			std::setfill(' ') << std::setw(10) << utility::toPrettyDiskString(writer.stream.tellp()) << '\t' <<
			std::setfill(' ') << std::setw(8) << (double)reader.stream.tellg()/reader.filesize*100 << "%" << ' ' <<
			timer.ElapsedString() << ' ' <<
			header.contigs[reader.front().body->CHROM].name << ":" << reader.front().body->POS+1 << "->" << reader.back().body->POS+1 << std::endl;
		}

		// Reset and update
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

		U64 total_uncompressed = 0; U64 total_compressed = 0;
		for(U32 i = 0; i < 10; ++i){
			total_uncompressed += this->stats_basic[i].cost_uncompressed;
			total_compressed   += this->stats_basic[i].cost_compressed;
		}

		std::cerr << utility::timestamp("PROGRESS") << "Wrote: " << utility::ToPrettyString(this->writer.n_variants_written) << " variants in " << utility::ToPrettyString(this->writer.n_blocks_written) << " blocks in " << timer.ElapsedString() << " to " << utility::toPrettyDiskString((U64)this->writer.stream.tellp()) << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << "BCF: " << utility::toPrettyDiskString(reader.filesize) << "\t" << utility::toPrettyDiskString(reader.b_data_read) << std::endl;
		std::cerr << utility::timestamp("PROGRESS") << "YON: " << utility::toPrettyDiskString(total_compressed) << '\t' << utility::toPrettyDiskString(total_uncompressed) << std::endl;
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

	meta_type meta(entry, this->block.header.minPosition);
	if(!this->parseBCFBody(meta, entry)){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode BCF body..." << std::endl;
		return false;
	}

	// GT encoding if available
	if(entry.hasGenotypes){
		meta.controller.gt_available = true;
		if(!this->encoder.Encode(entry, meta, this->block, this->permutator.manager->get())){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode GT..." << std::endl;
			return false;
		}
	} else {
		meta.controller.gt_available = false;
	}

	// Add meta
	this->block += meta;

	// Update number of entries in block
	++this->index_entry.n_variants;

	return true;
}

bool VariantImporter::parseBCFBody(meta_type& meta, bcf_entry_type& entry){
	for(U32 i = 0; i < entry.filterPointer; ++i){
		assert(entry.filterID[i].mapID != -1);
		this->block.AddFieldFILTER(this->header->filter_remap[entry.filterID[i].mapID]);
	}

	for(U32 i = 0; i < entry.infoPointer; ++i){
		assert(entry.infoID[i].mapID != -1);
		const U32 mapID = this->block.AddFieldINFO(this->header->info_remap[entry.infoID[i].mapID]);

		stream_container& target_container = this->block.info_containers[mapID];

		// Flags and integers
		// These are BCF value types
		U32 internal_pos = entry.infoID[i].l_offset;
		if(entry.infoID[i].primitive_type <= 3){
			for(U32 j = 0; j < entry.infoID[i].l_stride; ++j){
				target_container.Add(entry.getInteger(entry.infoID[i].primitive_type, internal_pos));
			}
		}
		// Floats
		else if(entry.infoID[i].primitive_type == bcf::BCF_FLOAT){
			for(U32 j = 0; j < entry.infoID[i].l_stride; ++j){
				target_container.Add(entry.getFloat(internal_pos));
			}
		}
		// Chars
		else if(entry.infoID[i].primitive_type == bcf::BCF_CHAR){
			target_container.AddCharacter(entry.getCharPointer(internal_pos), entry.infoID[i].l_stride);
			internal_pos += entry.infoID[i].l_stride;
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible in info: " << (int)entry.infoID[i].primitive_type << std::endl;
			exit(1);
		}

		++target_container;
		target_container.addStride(entry.infoID[i].l_stride);
	}

	for(U32 i = 0; i < entry.formatPointer; ++i){
		assert(entry.formatID[i].mapID != -1);

		//const U32 mapID = this->block.format_fields.setGet(this->header->format_remap[entry.formatID[i].mapID]);
		const U32 mapID = this->block.AddFieldFORMAT(this->header->format_remap[entry.formatID[i].mapID]);
		U32 internal_pos = entry.formatID[i].l_offset;

		// First value is always genotypes if there are any
		if(entry.hasGenotypes == true && i == 0)
			continue;

		// Hash INFO values
		stream_container& target_container = this->block.format_containers[mapID];

		// Flags and integers
		// These are BCF value types
		if(entry.formatID[i].primitive_type <= 3){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container.Add(entry.getInteger(entry.formatID[i].primitive_type, internal_pos));
			}
		}
		// Floats
		else if(entry.formatID[i].primitive_type == bcf::BCF_FLOAT){
			for(U32 s = 0; s < this->header->samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container.Add(entry.getFloat(internal_pos));
			}
		}
		// Chars
		else if(entry.formatID[i].primitive_type == bcf::BCF_CHAR){
			for(U32 s = 0; s < this->header->samples; ++s){
				//for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
				target_container.AddCharacter(entry.getCharPointer(internal_pos), entry.formatID[i].l_stride);
				internal_pos += entry.formatID[i].l_stride;
			}
		}
		// Illegal: parsing error
		else {
			std::cerr << "impossible: " << (int)entry.formatID[i].primitive_type << std::endl;
			std::cerr << utility::timestamp("LOG") << entry.formatID[i].mapID << '\t' << entry.formatID[i].l_stride << '\t' << (int)entry.formatID[i].primitive_type << '\t' << internal_pos << '/' << entry.l_data << std::endl;
			exit(1);
		}

		++target_container;
		target_container.addStride(entry.formatID[i].l_stride);
	}

	if(entry.filterPointer){
		// Hash FILTER pattern
		const U64 hash_filter_vector = entry.hashFilter();

		S32 mapID = this->block.getPatternsFILTER(hash_filter_vector);
		if(mapID == -1){
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.filterPointer; ++i)
				ret_pattern.push_back(this->header->filter_remap[entry.filterID[i].mapID]);

			mapID = this->block.filter_patterns.size();
			assert(mapID < 65536);
			this->block.addPatternFILTER(ret_pattern, hash_filter_vector);
		}
		meta.filter_pattern_id = mapID;
	}

	if(entry.infoPointer){
		// Hash INFO pattern
		const U64 hash_info_vector = entry.hashInfo();

		S32 mapID = this->block.getPatternsINFO(hash_info_vector);
		if(mapID == -1){
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.infoPointer; ++i)
				ret_pattern.push_back(this->header->info_remap[entry.infoID[i].mapID]);

			mapID = this->block.info_patterns.size();
			assert(mapID < 65536);
			this->block.addPatternINFO(ret_pattern, hash_info_vector);
		}
		meta.info_pattern_id = mapID;
	}

	if(entry.formatPointer){
		// Hash FORMAT pattern
		const U64 hash_format_vector = entry.hashFormat();

		S32 mapID = this->block.getPatternsFORMAT(hash_format_vector);
		if(mapID == -1){
			std::vector<U32> ret_pattern;
			for(U32 i = 0; i < entry.formatPointer; ++i)
				ret_pattern.push_back(this->header->format_remap[entry.formatID[i].mapID]);

			mapID = this->block.format_patterns.size();
			assert(mapID < 65536);
			this->block.addPatternFORMAT(ret_pattern, hash_format_vector);
		}
		meta.format_pattern_id = mapID;
	}

	// Return
	return true;
}

} /* namespace Tachyon */
