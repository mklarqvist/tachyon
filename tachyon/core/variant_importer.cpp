#include <fstream>
#include <regex>

#include "../algorithm/digital_digest.h"
#include "footer/footer.h"
#include "../containers/checksum_container.h"
#include "variant_importer.h"
#include "../algorithm/encryption/EncryptionDecorator.h"

namespace tachyon {

#define IMPORT_ASSERT 1

S32 reg2bin(S32 beg, S32 end){
	//--end;
	if (beg>>11 == end>>11) return ((1<<12)-1)/7 + (beg>>11);
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
	return 0;
}

S32 reg2bin2(S64 beg, S64 end, U64 ref_size, BYTE n_levels){
	std::cerr << "input: " << beg << "-" << end << std::endl;
	for(S32 i = n_levels; i >= 0; --i){
		std::cerr << i << ": " << ref_size / pow(4,i) << " " << S64(beg/(ref_size / pow(4,i))) << " == " << S64(end/(ref_size / pow(4,i))) << std::endl;
	}
	return(0);
}

VariantImporter::VariantImporter(std::string inputFile,
		                         std::string outputPrefix,
                                   const U32 checkpoint_n_snps,
                                const double checkpoint_bases) :
	GT_available_(false),
	permute(true),
	encrypt(false),
	checkpoint_n_snps(checkpoint_n_snps),
	checkpoint_bases(checkpoint_bases),
	info_end_key_(-1),
	info_svlen_key_(-1),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	writer(nullptr),
	header(nullptr)
{
}

VariantImporter::~VariantImporter(){
	delete this->writer;
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

	encryption::EncryptionDecorator encryptionManager;
	encryption::Keychain keychain;

	this->header = &reader.header;

	// Spawn RLE controller and update PPA controller
	this->encoder.setSamples(this->header->samples);
	this->block.ppa_manager.setSamples(this->header->samples);
	this->permutator.manager = &this->block.ppa_manager;
	this->permutator.setSamples(this->header->samples);

	if(this->outputPrefix.size() == 0) this->writer = new writer_stream_type;
	else this->writer = new writer_file_type;

	if(!this->writer->open(this->outputPrefix)){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	// Writer MAGIC
	this->writer->stream->write(&tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH);
	// Convert VCF header to Tachyon heeader
	core::VariantHeader header(*this->header);
	header.literals += "\n##tachyon_importVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();
	header.literals += "\n##tachyon_importCommand=import -i " + this->inputFile + " -o " + this->outputPrefix + " -c " + std::to_string(this->checkpoint_n_snps) + " -C " + std::to_string(this->checkpoint_bases);
	if(this->encrypt) header.literals += " -k";
	if(this->permute) header.literals += " -P";
	else header.literals += " -p";
	header.header_magic.l_literals = header.literals.size();

	// Convert header to byte stream, compress, and write to file
	containers::DataContainer header_data;
	header_data.resize(65536 + header.literals.size()*2);
	header_data.buffer_data_uncompressed << header;
	this->compression_manager.zstd_codec.compress(header_data);
	*this->writer->stream << header_data.header; // write header
	*this->writer->stream << header_data.buffer_data;

	// Search for GT field in the header
	this->GT_available_ = header.has_format_field("GT");
	for(U32 i = 0; i < this->header->format_map.size(); ++i){
		if(this->header->format_map[i].ID == "GT"){
			reader.map_gt_id = this->header->format_map[i].IDX;
		}
	}

	// Search for END field in the header
	for(U32 i = 0; i < this->header->info_map.size(); ++i){
		if(this->header->info_map[i].ID == "END"){
			this->info_end_key_ = this->header->info_map[i].IDX;
			std::cerr << "Found END at: " << this->header->info_map[i].IDX << std::endl;
		}
	}

	// Search for END field in the header
	for(U32 i = 0; i < this->header->info_map.size(); ++i){
		if(this->header->info_map[i].ID == "SVLEN"){
			this->info_svlen_key_ = this->header->info_map[i].IDX;
			std::cerr << "Found SVLEN at: " << this->header->info_map[i].IDX << std::endl;
		}
	}

	this->block.header.controller.hasGT = this->GT_available_;

	// Resize containers
	const U32 resize_to = this->checkpoint_n_snps * sizeof(U32) * 2; // small initial allocation
	this->block.resize(resize_to);

	// Digest controller
	algorithm::VariantDigitalDigestManager checksums(25, this->header->info_map.size(), this->header->format_map.size());

	// Start import
	U32 previousFirst    = 0;
	U32 previousLast     = 0;
	S32 previousContigID = -1;

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

		// Update head meta
		this->block.header.controller.hasGT = this->GT_available_;
		this->block.header.n_variants       = reader.size();
		this->block.finalize();

		// Perform compression using standard parameters
		if(!this->compression_manager.compress(this->block)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Checksum have to come before encryption
		checksums += this->block;

		// Encryption
		if(this->encrypt){
			this->block.header.controller.anyEncrypted = true;
			if(!encryptionManager.encrypt(this->block, keychain, YON_ENCRYPTION_AES_256_GCM)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to encrypt..." << std::endl;
			}
		}

		this->index_entry.byte_offset = this->writer->stream->tellp();
		this->block.write(*this->writer->stream, this->stats_basic, this->stats_info, this->stats_format);

		// Compress and write footer
		this->block.footer_support.buffer_data_uncompressed << this->block.footer;
		this->compression_manager.zstd_codec.compress(this->block.footer_support);
		const U64 start_footer_pos = this->writer->stream->tellp();
		const U32 footer_uLength   = this->block.footer_support.header.data_header.uLength;
		const U32 footer_cLength   = this->block.footer_support.header.data_header.cLength;
		const U32 footer_crc       = this->block.footer_support.header.data_header.crc;
		this->writer->stream->write(reinterpret_cast<const char*>(&footer_uLength), sizeof(U32));
		this->writer->stream->write(reinterpret_cast<const char*>(&footer_cLength), sizeof(U32));
		this->writer->stream->write(reinterpret_cast<const char*>(&footer_crc),     sizeof(U32));
		*this->writer->stream << this->block.footer_support.buffer_data;

		stats_basic[0].cost_uncompressed += (U64)this->writer->stream->tellp() - start_footer_pos;

		// Write EOB
		this->writer->stream->write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));

		this->index_entry.blockID         = this->block.header.blockID;
		this->index_entry.byte_offset_end = this->writer->stream->tellp();
		this->index_entry.contigID        = reader.front().body->CHROM;
		this->index_entry.minPosition     = reader.front().body->POS;
		this->index_entry.maxPosition     = reader.back().body->POS;
		this->index_entry.n_variants      = reader.size();
		this->writer->index += this->index_entry;

		this->index_entry.reset();
		++this->writer->n_blocks_written;
		this->writer->n_variants_written += reader.size();

		if(!SILENT){
			std::cerr << utility::timestamp("PROGRESS") <<
			std::setfill(' ') << std::setw(10) << this->writer->n_variants_written << ' ' <<
			std::setfill(' ') << std::setw(10) << utility::toPrettyDiskString(this->writer->stream->tellp()) << '\t' <<
			std::setfill(' ') << std::setw(8)  << (double)reader.stream.tellg()/reader.filesize*100 << "%" << ' ' <<
			timer.ElapsedString() << ' ' <<
			header.contigs[reader.front().body->CHROM].name << ":" << reader.front().body->POS+1 << "->" << reader.back().body->POS+1 << std::endl;
		}

		// Reset and update
		this->block.clear();
		this->permutator.reset();
		this->writer->stream->flush();
		previousContigID = reader.front().body->CHROM;
		previousFirst    = reader.front().body->POS;
		previousLast     = reader.back().body->POS;
		this->index_entry.reset();
	}
	//t Done importing
	this->writer->stream->flush();

	core::Footer footer;
	footer.offset_end_of_data = this->writer->stream->tellp();
	footer.n_blocks           = this->writer->n_blocks_written;
	footer.n_variants         = this->writer->n_variants_written;
	assert(footer.n_blocks == this->writer->index.size());


	this->writer->writeIndex(); // Write index
	checksums.finalize();       // Finalize SHA-512 digests
	*this->writer->stream << checksums;
	*this->writer->stream << footer;
	this->writer->stream->flush();

		std::vector<std::string> usage_statistics_names = {
			"FooterHeader","GT-PPA","MetaContig","MetaPositions","MetaRefAlt","MetaController","MetaQuality","MetaNames",
			"MetaAlleles","MetaInfoMaps","MetaFormatMaps","MetaFilterMaps","GT-Support",
			"GT-RLE8","GT-RLE16","GT-RLE32","GT-RLE64",
			"GT-Simple8","GT-Simple16","GT-Simple32","GT-Simple64","INFO-ALL","FORMAT-ALL"};

		U64 total_uncompressed = 0; U64 total_compressed = 0;
		for(U32 i = 0; i < usage_statistics_names.size(); ++i){
			total_uncompressed += this->stats_basic[i].cost_uncompressed;
			total_compressed   += this->stats_basic[i].cost_compressed;
		}

		if(this->outputPrefix.size()){
			std::ofstream writer_stats;
			writer_file_type* wstats = reinterpret_cast<writer_file_type*>(this->writer);
			writer_stats.open(wstats->basePath + wstats->baseName + "_yon_stats.txt", std::ios::out);

			if(!SILENT)
				std::cerr << utility::timestamp("LOG") << "Writing statistics to: " << (wstats->basePath + wstats->baseName) << "_yon_stats.txt" << std::endl;

			if(writer_stats.good()){
				for(U32 i = 0; i < usage_statistics_names.size(); ++i) writer_stats << usage_statistics_names[i] << '\t' << this->stats_basic[i] << std::endl;
				for(U32 i = 0; i < header.header_magic.n_info_values; ++i) writer_stats << "INFO_" << header.info_fields[i].ID << '\t' << this->stats_info[i] << std::endl;
				for(U32 i = 0; i < header.header_magic.n_format_values; ++i) writer_stats << "FORMAT_" << header.format_fields[i].ID << '\t' << this->stats_format[i] << std::endl;

				writer_stats << "BCF\t" << reader.filesize << "\t" << reader.b_data_read << '\t' << (float)reader.b_data_read/reader.filesize << std::endl;
				writer_stats << "YON\t" << this->writer->stream->tellp() << "\t" << total_uncompressed << '\t' << (float)reader.b_data_read/this->writer->stream->tellp() << std::endl;
				writer_stats.close();
			} else {
				std::cerr << utility::timestamp("ERROR", "SUPPORT")  << "Failed to open: " << (wstats->basePath + wstats->baseName + "_yon_stats.txt") << "... Continuing..." << std::endl;
			}

			const algorithm::GenotypeEncoderStatistics& gt_stats = this->encoder.getUsageStats();
			const U64 n_total_gt = gt_stats.getTotal();
			if(!SILENT){
				std::cout << "GT-RLE-8\t"   << gt_stats.rle_counts[0] << '\t' << (float)gt_stats.rle_counts[0]/n_total_gt << std::endl;
				std::cout << "GT-RLE-16\t"  << gt_stats.rle_counts[1] << '\t' << (float)gt_stats.rle_counts[1]/n_total_gt << std::endl;
				std::cout << "GT-RLE-32\t"  << gt_stats.rle_counts[2] << '\t' << (float)gt_stats.rle_counts[2]/n_total_gt << std::endl;
				std::cout << "GT-RLE-64\t"  << gt_stats.rle_counts[3] << '\t' << (float)gt_stats.rle_counts[3]/n_total_gt << std::endl;
				std::cout << "GT-RLES-8\t"  << gt_stats.rle_simple_counts[0] << '\t' << (float)gt_stats.rle_simple_counts[0]/n_total_gt << std::endl;
				std::cout << "GT-RLES-16\t" << gt_stats.rle_simple_counts[1] << '\t' << (float)gt_stats.rle_simple_counts[1]/n_total_gt << std::endl;
				std::cout << "GT-RLES-32\t" << gt_stats.rle_simple_counts[2] << '\t' << (float)gt_stats.rle_simple_counts[2]/n_total_gt << std::endl;
				std::cout << "GT-RLES-64\t" << gt_stats.rle_simple_counts[3] << '\t' << (float)gt_stats.rle_simple_counts[3]/n_total_gt << std::endl;
				std::cout << "GT-DIPLOID-BCF-8\t"  << gt_stats.diploid_bcf_counts[0] << '\t' << (float)gt_stats.diploid_bcf_counts[0]/n_total_gt << std::endl;
				std::cout << "GT-DIPLOID-BCF-16\t" << gt_stats.diploid_bcf_counts[1] << '\t' << (float)gt_stats.diploid_bcf_counts[1]/n_total_gt << std::endl;
				std::cout << "GT-DIPLOID-BCF-32\t" << gt_stats.diploid_bcf_counts[2] << '\t' << (float)gt_stats.diploid_bcf_counts[2]/n_total_gt << std::endl;
				std::cout << "GT-BCF-8\t"  << gt_stats.bcf_counts[0] << '\t' << (float)gt_stats.bcf_counts[0]/n_total_gt << std::endl;
				std::cout << "GT-BCF-16\t" << gt_stats.bcf_counts[1] << '\t' << (float)gt_stats.bcf_counts[1]/n_total_gt << std::endl;
				std::cout << "GT-BCF-32\t" << gt_stats.bcf_counts[2] << '\t' << (float)gt_stats.bcf_counts[2]/n_total_gt << std::endl;
				std::cerr << utility::timestamp("PROGRESS") << "Wrote: " << utility::ToPrettyString(this->writer->n_variants_written) << " variants in " << utility::ToPrettyString(this->writer->n_blocks_written) << " blocks in " << timer.ElapsedString() << " to " << utility::toPrettyDiskString((U64)this->writer->stream->tellp()) << std::endl;
			}
		}

	if(this->encrypt){
		if(this->outputPrefix.size()){
			std::ofstream writer_keychain;
			writer_file_type* wstats = reinterpret_cast<writer_file_type*>(this->writer);
			writer_keychain.open(wstats->basePath + wstats->baseName + ".kyon", std::ios::out);
			if(!SILENT)
				std::cerr << utility::timestamp("LOG") << "Writing encryption keychain to: " << (wstats->basePath + wstats->baseName) << ".kyon" << std::endl;

			if(writer_keychain.good()){
				writer_keychain.write(constants::FILE_HEADER.data(), constants::FILE_HEADER_LENGTH);
				writer_keychain.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MAJOR),   sizeof(S32));
				writer_keychain.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MINOR),   sizeof(S32));
				writer_keychain.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_RELEASE), sizeof(S32));

				//writer_keychain.write(keybuffer.data(), keybuffer.size());
				writer_keychain << keychain;
				writer_keychain.flush();
			}
			const U32 keychain_size = writer_keychain.tellp();
			writer_keychain.close();
			if(!SILENT)
				std::cerr << utility::timestamp("LOG") << "Wrote keychain with " << utility::ToPrettyString(keychain.size()) << " keys to " << utility::toPrettyDiskString(keychain_size) << "..." << std::endl;
		}
	}

	// temp
	for(U32 i = 0; i < writer->index.size(); ++i){
		std::cerr << i << "\t";
		writer->index.getIndex().at(i).print(std::cerr);
		std::cerr << std::endl;
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

	// Has END position?
	S32 index_bin = -1;

	/*
	// test
	//std::cerr << this->header->contigs[entry.body->CHROM].bp_length << std::endl;
	const U64 contig_length = this->header->contigs[entry.body->CHROM].bp_length;
	const double raw_section_range = (double)contig_length/pow(4,8);
	const U64 used_contig_length = ceil(ceil(raw_section_range) / ((double)contig_length/pow(4,8))*contig_length);
	//std::cerr << raw_section_range << "->" << used_contig_length << " sub " << ((double)contig_length/pow(4,8))*contig_length << std::endl;

	//for(S32 i = 7; i > 0; --i){
	//	std::cerr << i << ' ' << pow(4, i) << '\t' << (float)used_contig_length/pow(4,i) << std::endl;
	//}


	if(this->info_end_key_ != -1){
		// Linear search: this is not optimal but probably still faster
		// than generating a new hash table for each record
		for(U32 i = 0; i < entry.infoPointer; ++i){
			if(entry.infoID[i].mapID == this->info_end_key_){
				const U32 end = entry.getInteger(entry.infoID[i].primitive_type, entry.infoID[i].l_offset);
				std::cerr << "Found END at " << i << ".  END=" << end << " POS=" << entry.body->POS+1 << " BIN=" << reg2bin(entry.body->POS, end)  << std::endl;
				reg2bin2(entry.body->POS, end, used_contig_length, 7);

				index_bin = reg2bin(entry.body->POS, end);
				break;
			}
		}
	}

	if(index_bin == -1){
		S32 longest = -1;
		for(U32 i = 0; i < meta.n_alleles; ++i){
			// Regex pattern ^[ATGC]{1,}$
			std::regex txt_regex("^[ATGC]{1,}$");
			if(std::regex_match(meta.alleles[i].toString(), txt_regex)){
				//std::cerr << "regex match: " << meta.alleles[i].toString() << std::endl;
				if(meta.alleles[i].l_allele > longest) longest = meta.alleles[i].l_allele;
			} else {
				std::cerr << "POS=" << entry.body->POS+1 << " no regex match: " << meta.alleles[i].toString() << std::endl;
			}
			// alt = <> or
			//meta.alleles[i].l_allele
		}

		if(longest){
			//if(longest > 3) std::cerr << "longest: " << longest << std::endl;
			index_bin = reg2bin(entry.body->POS, entry.body->POS + longest);
		} else // fallback if all others fail
			index_bin = reg2bin(entry.body->POS, entry.body->POS);
	}
	*/
	index_bin = reg2bin(entry.body->POS, entry.body->POS);
	if(index_bin > this->index_entry.maxBin) this->index_entry.maxBin = index_bin;
	if(index_bin < this->index_entry.minBin) this->index_entry.minBin = index_bin;

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
