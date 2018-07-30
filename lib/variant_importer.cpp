#include <fstream>
#include <regex>

#include "variant_importer.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "containers/checksum_container.h"
#include "core/footer/footer.h"

namespace tachyon {

VariantImporter::VariantImporter(const settings_type& settings) :
		settings_(settings),
		GT_available_(false),
		writer(nullptr),
		header(nullptr)
{

}

VariantImporter::~VariantImporter(){
	delete this->writer;
}

void VariantImporter::clear(){
	this->vcf_container_.clear();
}

bool VariantImporter::Build(){
	std::ifstream temp(this->settings_.input_file, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << utility::timestamp("ERROR", "IMPORT")  << "Failed to open file (" << this->settings_.input_file << ")..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if((BYTE)tempData[0] == io::constants::GZIP_ID1 && (BYTE)tempData[1] == io::constants::GZIP_ID2){
		// if(!this->BuildBCF()){
		if(!this->BuildVCF()){
			std::cerr << utility::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << utility::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool VariantImporter::BuildVCF(void){
	// Retrieve a unique VcfReader.
	this->vcf_reader_ = io::VcfReader::FromFile(this->settings_.input_file);
	if(this->vcf_reader_ == nullptr){
		return false;
	}

	// Remap the global IDX fields in Vcf to the appropriate incremental order.
	// This is useful in the situations when fields have been removed or added
	// to the Vcf header section without reformatting the file.
	for(U32 i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i)
		this->contig_reorder_map_[this->vcf_reader_->vcf_header_.contigs_[i].idx] = i;

	for(U32 i = 0; i < this->vcf_reader_->vcf_header_.info_fields_.size(); ++i)
		this->info_reorder_map_[this->vcf_reader_->vcf_header_.info_fields_[i].idx] = i;

	for(U32 i = 0; i < this->vcf_reader_->vcf_header_.format_fields_.size(); ++i)
		this->format_reorder_map_[this->vcf_reader_->vcf_header_.format_fields_[i].idx] = i;

	for(U32 i = 0; i < this->vcf_reader_->vcf_header_.filter_fields_.size(); ++i)
		this->filter_reorder_map_[this->vcf_reader_->vcf_header_.filter_fields_[i].idx] = i;

	// Predicate of a search for "GT" FORMAT field in the Vcf header.
	this->GT_available_ = (this->vcf_reader_->vcf_header_.GetFormat("GT") != nullptr);

	// Predicate of a search for "END" INFO field in the Vcf header.
	io::VcfInfo* vcf_info_end = this->vcf_reader_->vcf_header_.GetInfo("END");
	if(vcf_info_end != nullptr)
		this->settings_.info_end_key = vcf_info_end->idx;

	// Allocate a new writer.
	if(this->settings_.output_prefix.size() == 0 ||
	   (this->settings_.output_prefix.size() == 1 && this->settings_.output_prefix == "-"))
	{
		this->writer = new writer_stream_type;
	}
	else this->writer = new writer_file_type;

	// Setup the checksums container.
	algorithm::VariantDigestManager checksums(YON_BLK_N_STATIC + 1, // Add one for global checksum.
			this->vcf_reader_->vcf_header_.info_fields_.size() + 1,
			this->vcf_reader_->vcf_header_.format_fields_.size() + 1);

	// The index needs to know how many contigs that's described in the
	// Vcf header and their lenghts. This information is needed to construct
	// the linear and quad-tree index most appropriate for the data.
	this->writer->index.Add(this->vcf_reader_->vcf_header_.contigs_);

	// Resize containers
	const U32 resize_to = this->settings_.checkpoint_n_snps * sizeof(U32) * 2; // small initial allocation
	this->block.resize(resize_to);

	// Todo:
	this->permutator.SetSamples(this->vcf_reader_->vcf_header_.GetNumberSamples());

	uint64_t b_total_genotypes = 0;
	uint64_t b_total_write = 0;

	// Iterate over all available variants in the file or until encountering
	// an error.
	while(true){
		this->block.footer.AllocateHeaders(
				this->vcf_reader_->vcf_header_.info_fields_.size(),
				this->vcf_reader_->vcf_header_.format_fields_.size(),
				this->vcf_reader_->vcf_header_.filter_fields_.size());

		// Retrieve bcf1_t records using htslib and lazy evaluate them. Stop
		// after retrieving a set number of variants or if the interval between
		// the smallest and largest variant exceeds some distance in base pairs.
		if(this->vcf_container_.getVariants(this->settings_.checkpoint_n_snps,
		                                    this->settings_.checkpoint_bases,
		                                    this->vcf_reader_) == false){
			break;
		}

		std::cerr << this->vcf_container_.sizeWithoutCarryOver() << ": " <<
		             this->vcf_container_.front()->pos + 1 << "->" <<
					 this->vcf_container_.back()->pos + 1 << std::endl;

		// Todo
		if(this->permutator.Build(this->vcf_container_) == false) return false;
		if(this->AddRecords(this->vcf_container_) == false) return false;

		this->block.header.controller.hasGT  = this->GT_available_;
		this->block.UpdateContainers();
		this->block.Finalize();

		// Perform compression using standard parameters
		if(!this->compression_manager.compress(this->block, this->settings_.compression_level, 6)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Update checksums container with the available data.
		checksums += this->block;

		// Todo: this->index += this->vcf_reader_->variant_reader

		uint64_t n_sum = 0;
		for(U32 i = YON_BLK_GT_INT8; i < YON_BLK_GT_S_INT64; ++i){
			//std::cerr << this->block.base_containers[i].buffer_data_uncompressed.size() << "->" << this->block.base_containers[i].buffer_data.size() << "=" << (float)this->block.base_containers[i].buffer_data_uncompressed.size()/this->block.base_containers[i].buffer_data.size() << "\t";
			n_sum += this->block.base_containers[i].buffer_data_uncompressed.size();
			b_total_genotypes += this->block.base_containers[i].buffer_data.size();
		}
		//std::cerr << this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data_uncompressed.size() << "->" << this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data.size() << "=" << (float)this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data_uncompressed.size()/this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data.size() << std::endl;
		n_sum += this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data_uncompressed.size();
		b_total_genotypes += this->block.base_containers[YON_BLK_GT_S_INT64].buffer_data.size();
		b_total_genotypes += this->block.base_containers[YON_BLK_GT_SUPPORT].buffer_data.size();
		b_total_write += this->block.DetermineCompressedSize();
		std::cerr << "Total: " << b_total_genotypes << " local: " << n_sum << std::endl;
		std::cerr << "PPA size: " << this->block.ppa_manager.getObjectSize() << std::endl;
		std::cerr << "Written: " << this->block.DetermineCompressedSize() << " -> " << b_total_write << std::endl;

		// Clear current data.
		this->clear();
		this->block.clear();
		this->index_entry.reset();
	}
	std::cerr << "Total write: " << b_total_write << std::endl;

	// Finalize checksum container.
	checksums.Finalize();

	// All done
	return(true);
}

bool VariantImporter::BuildBCF(void){
	bcf_reader_type bcf_reader;
	if(!bcf_reader.open(this->settings_.input_file)){
		std::cerr << utility::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	encryption::EncryptionDecorator encryption_manager;
	encryption::Keychain<> keychain;

	this->header = &bcf_reader.header;

	// Spawn RLE controller and update PPA controller
	this->encoder.setSamples(this->header->n_samples);
	this->block.ppa_manager.setSamples(this->header->n_samples);
	this->permutator.manager = &this->block.ppa_manager;
	this->permutator.SetSamples(this->header->n_samples);

	if(this->settings_.output_prefix.size() == 0) this->writer = new writer_stream_type;
	else this->writer = new writer_file_type;

	if(!this->writer->open(this->settings_.output_prefix)){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	//index::VariantIndex idx;
	for(U32 i = 0; i < this->header->contigs.size(); ++i){
		const U64 contig_length = this->header->contigs[i].bp_length;
		BYTE n_levels = 7;
		U64 bins_lowest = pow(4,n_levels);
		double used = ( bins_lowest - (contig_length % bins_lowest) ) + contig_length;

		if(used / bins_lowest < 2500){
			for(S32 i = n_levels; i != 0; --i){
				if(used/pow(4,i) > 2500){
					n_levels = i;
					break;
				}
			}
		}

		this->writer->index.index_.add(i, contig_length, n_levels);
		//std::cerr << "contig: " << this->header->contigs[i].name << "(" << i << ")" << " -> " << contig_length << " levels: " << (int)n_levels << std::endl;
		//std::cerr << "idx size:" << idx.size() << " at " << this->writer->index.variant_index_[i].size() << std::endl;
		//std::cerr << i << "->" << this->header->contigs[i].name << ":" << contig_length << " up to " << (U64)used << " width (bp) lowest level: " << used/pow(4,n_levels) << "@level: " << (int)n_levels << std::endl;
	}

	// Writer MAGIC
	this->writer->stream->write(&tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH);
	// Convert VCF header to Tachyon heeader
	core::VariantHeader header(*this->header);
	header.literals += "\n##tachyon_importVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
	                +  SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();
	header.literals += "\n" + this->settings_.GetInterpretedString();

	if(this->settings_.encrypt_data) header.literals += " -k";
	if(this->settings_.permute_genotypes) header.literals += " -P";
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
			bcf_reader.map_gt_id = this->header->format_map[i].IDX;
		}
	}

	// Search for END field in the header
	for(U32 i = 0; i < this->header->info_map.size(); ++i){
		if(this->header->info_map[i].ID == "END"){
			this->settings_.info_end_key = this->header->info_map[i].IDX;
			//std::cerr << "Found END at: " << this->header->info_map[i].IDX << std::endl;
		}
	}

	// Search for END field in the header
	for(U32 i = 0; i < this->header->info_map.size(); ++i){
		if(this->header->info_map[i].ID == "SVLEN"){
			this->settings_.info_svlen_key = this->header->info_map[i].IDX;
			//std::cerr << "Found SVLEN at: " << this->header->info_map[i].IDX << std::endl;
		}
	}

	// Set flag if genotypes are available
	this->block.header.controller.hasGT = this->GT_available_;

	// Resize containers
	const U32 resize_to = this->settings_.checkpoint_n_snps * sizeof(U32) * 2; // small initial allocation
	this->block.resize(resize_to);

	// Digest controller
	algorithm::VariantDigestManager checksums(25, this->header->info_map.size(), this->header->format_map.size());

	// Start import
	U32 previous_first     = 0;
	U32 previous_last      = 0;
	S32 previous_contig_ID = -1;

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

	bcf_reader.setFilterInvariantSites(this->settings_.drop_invariant_sites);
	while(true){
		// Retrieve BCF records
		if(!bcf_reader.getVariants(this->settings_.checkpoint_n_snps, this->settings_.checkpoint_bases)){
			break;
		}

		// Debug assertion
#if IMPORT_ASSERT == 1
		if(bcf_reader.front().body->CHROM == previous_contig_ID){
			if(!(bcf_reader.front().body->POS >= previous_first && bcf_reader.front().body->POS >= previous_last)){
				std::cerr << utility::timestamp("ERROR","IMPORT") << bcf_reader.front().body->POS << '/' << previous_first << '/' << previous_last << std::endl;
				std::cerr << bcf_reader[bcf_reader.n_entries].body->POS << std::endl;
				exit(1);
			}
		}
#endif
		this->block.header.blockID     = this->writer->n_blocks_written;
		this->block.header.contigID    = bcf_reader.front().body->CHROM;
		this->block.header.minPosition = bcf_reader.front().body->POS;
		this->block.header.maxPosition = bcf_reader.back().body->POS; // correct only for SNVs
		this->block.header.controller.hasGT         = this->GT_available_;
		this->block.header.controller.hasGTPermuted = this->settings_.permute_genotypes;
		// if there is 0 or 1 samples then GT data is never permuted
		if(header.getSampleNumber() <= 1)
			this->block.header.controller.hasGTPermuted = false;

		// Permute GT if GT is available and the appropriate flag is triggered
		if(this->block.header.controller.hasGT && this->block.header.controller.hasGTPermuted){
			if(!this->permutator.Build(bcf_reader)){
				std::cerr << utility::timestamp("ERROR","PERMUTE") << "Failed to complete..." << std::endl;
				return false;
			}
		}

		//\////////////////////////////////////////////////
		// Start new
		// Perform parsing of BCF entries in memory
		// Split out RLE compression and importing other INFO/FORMAT
		meta_type* meta_entries = static_cast<meta_type*>(::operator new[](bcf_reader.size() * sizeof(meta_type)));

		// Load meta data
		for(U32 i = 0; i < bcf_reader.size(); ++i){
			new( meta_entries + i ) meta_type( bcf_reader[i], this->block.header.minPosition );
			if(!this->addSite(meta_entries[i], bcf_reader[i])){
				std::cerr << utility::timestamp("ERROR","IMPORT") << "Failed to add BCF entry..." << std::endl;
				return false;
			}
		}
		// Add genotypes in parallel
		this->addGenotypes(bcf_reader, meta_entries);
		// Overload
		for(U32 i = 0; i < bcf_reader.size(); ++i) this->block += meta_entries[i];

		// Clean up
		for(std::size_t i = 0; i < bcf_reader.size(); ++i) (meta_entries + i)->~MetaEntry();
		::operator delete[](static_cast<void*>(meta_entries));
		//\////////////////////////////////////////////////

		// Update head meta
		this->block.header.controller.hasGT = this->GT_available_;
		this->block.header.n_variants       = bcf_reader.size();
		this->block.finalize();

		// Perform compression using standard parameters
		if(!this->compression_manager.compress(this->block, this->settings_.compression_level, 6)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Checksum have to come before encryption
		checksums += this->block;

		// Encryption
		if(this->settings_.encrypt_data){
			this->block.header.controller.anyEncrypted = true;
			if(!encryption_manager.encrypt(this->block, keychain, YON_ENCRYPTION_AES_256_GCM)){
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

		// Update index
		this->index_entry.blockID         = this->block.header.blockID;
		this->index_entry.byte_offset_end = this->writer->stream->tellp();
		this->index_entry.contigID        = bcf_reader.front().body->CHROM;
		this->index_entry.minPosition     = bcf_reader.front().body->POS;
		//this->index_entry.maxPosition     = bcf_reader.back().body->POS;
		this->index_entry.n_variants      = bcf_reader.size();
		this->writer->index.index_.linear_at(index_entry.contigID) += this->index_entry; // Todo: beautify
		this->index_entry.reset();

		++this->writer->n_blocks_written;
		this->writer->n_variants_written += bcf_reader.size();
		++this->writer->index.number_blocks;

		if(!SILENT){
			std::cerr << utility::timestamp("PROGRESS") <<
			std::setfill(' ') << std::setw(10) << this->writer->n_variants_written << ' ' <<
			std::setfill(' ') << std::setw(10) << utility::toPrettyDiskString(this->writer->stream->tellp()) << '\t' <<
			std::setfill(' ') << std::setw(8)  << (double)bcf_reader.stream.tellg()/bcf_reader.filesize*100 << "%" << ' ' <<
			timer.ElapsedString() << ' ' <<
			header.contigs[bcf_reader.front().body->CHROM].name << ":" << bcf_reader.front().body->POS+1 << "->" << bcf_reader.back().body->POS+1 << std::endl;
		}

		// Reset and update
		this->block.clear();
		this->permutator.reset();
		this->writer->stream->flush();
		previous_contig_ID = bcf_reader.front().body->CHROM;
		previous_first    = bcf_reader.front().body->POS;
		previous_last     = bcf_reader.back().body->POS;
		this->index_entry.reset();
	}
	// Done importing
	this->writer->stream->flush();

	core::Footer footer;
	footer.offset_end_of_data = this->writer->stream->tellp();
	footer.n_blocks           = this->writer->n_blocks_written;
	footer.n_variants         = this->writer->n_variants_written;
	assert(footer.n_blocks == this->writer->index.size());

	U64 last_pos = this->writer->stream->tellp();
	this->writer->writeIndex(); // Write index
	std::cerr << utility::timestamp("PROGRESS") << "Index size: " << utility::toPrettyDiskString((U64)this->writer->stream->tellp() - last_pos) << "..." << std::endl;
	last_pos = this->writer->stream->tellp();
	checksums.finalize();       // Finalize SHA-512 digests
	*this->writer->stream << checksums;
	std::cerr << utility::timestamp("PROGRESS") << "Checksum size: " << utility::toPrettyDiskString((U64)this->writer->stream->tellp() - last_pos) << "..." << std::endl;
	last_pos = this->writer->stream->tellp();
	*this->writer->stream << footer;
	std::cerr << utility::timestamp("PROGRESS") << "Footer size: " << utility::toPrettyDiskString((U64)this->writer->stream->tellp() - last_pos) << "..." << std::endl;

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

	// If we are writing to a file
	if(this->settings_.output_prefix.size()){
		std::ofstream writer_stats;
		writer_file_type* wstats = reinterpret_cast<writer_file_type*>(this->writer);
		writer_stats.open(wstats->basePath + wstats->baseName + "_yon_stats.txt", std::ios::out);

		if(!SILENT)
			std::cerr << utility::timestamp("LOG") << "Writing statistics to: " << (wstats->basePath + wstats->baseName) << "_yon_stats.txt" << std::endl;

		if(writer_stats.good()){
			for(U32 i = 0; i < usage_statistics_names.size(); ++i)       writer_stats << usage_statistics_names[i] << '\t' << this->stats_basic[i] << std::endl;
			for(U32 i = 0; i < header.header_magic.n_info_values; ++i)   writer_stats << "INFO_" << header.info_fields[i].ID << '\t' << this->stats_info[i] << std::endl;
			for(U32 i = 0; i < header.header_magic.n_format_values; ++i) writer_stats << "FORMAT_" << header.format_fields[i].ID << '\t' << this->stats_format[i] << std::endl;

			writer_stats << "BCF\t" << bcf_reader.filesize << "\t" << bcf_reader.b_data_read << '\t' << (float)bcf_reader.b_data_read/bcf_reader.filesize << std::endl;
			writer_stats << "YON\t" << this->writer->stream->tellp() << "\t" << total_uncompressed << '\t' << (float)bcf_reader.b_data_read/this->writer->stream->tellp() << std::endl;
			writer_stats.close();
		} else {
			std::cerr << utility::timestamp("ERROR", "SUPPORT")  << "Failed to open: " << (wstats->basePath + wstats->baseName + "_yon_stats.txt") << "... Continuing..." << std::endl;
		}
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

	// Write keychain
	if(this->settings_.encrypt_data){
		if(this->settings_.output_prefix.size()){
			std::ofstream writer_keychain;
			writer_file_type* wstats = reinterpret_cast<writer_file_type*>(this->writer);
			writer_keychain.open(wstats->basePath + wstats->baseName + ".kyon", std::ios::out);
			if(!SILENT)
				std::cerr << utility::timestamp("LOG") << "Writing encryption keychain to: " << (wstats->basePath + wstats->baseName) << ".kyon" << std::endl;

			if(writer_keychain.good()){
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
	// bins in contig
	/*
	for(U32 i = 0; i < this->writer->index.variant_index_.at(19).size(); ++i){
		std::cout << i << '\t' << this->writer->index.variant_index_.at(19).at(i).blockID << '\t' << this->writer->index.variant_index_.at(19).at(i).n_variants_ << '\t';
		// hits in bin
		for(U32 j = 0; j < this->writer->index.variant_index_.at(19).at(i).size(); ++j)
			std::cout << ',' << this->writer->index.variant_index_.at(19).at(i).at(j);
		std::cout << std::endl;
	}
	*/

	// All done
	return(true);
}

bool VariantImporter::addGenotypes(bcf_reader_type& bcf_reader, meta_type* meta_entries){
	/*
	for(U32 i = 0; i < bcf_reader.size(); ++i){
		if(bcf_reader[i].hasGenotypes){
			meta_entries[i].controller.gt_available = true;

			if(!this->encoder.Encode(bcf_reader[i], meta_entries[i], this->block, this->permutator.manager->get())){
				std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode GT..." << std::endl;
				return false;
			}
		} else {
			meta_entries[i].controller.gt_available = false;
		}
	}
	*/
	this->encoder.EncodeParallel(bcf_reader, meta_entries, this->block, this->permutator.manager->get(), this->settings_.n_threads);

	return true;
}

bool VariantImporter::AddRecords(const vcf_container_type& container){
	meta_type* meta_entries_new = static_cast<meta_type*>(::operator new[](container.sizeWithoutCarryOver() * sizeof(meta_type)));
	for(U32 i = 0; i < container.sizeWithoutCarryOver(); ++i){
		// Transmute a bcf record into a meta structure
		new( meta_entries_new + i ) meta_type( this->vcf_container_[i], this->block.header.minPosition );
		// Add the record data
		if(this->AddRecord(container, i, meta_entries_new[i]) == false)
			return false;
	}

	// Add genotypes in parallel
	// Todo: this->addGenotypes(bcf_reader, meta_entries);
	this->encoder.setSamples(this->vcf_reader_->vcf_header_.GetNumberSamples());
	this->encoder.EncodeNew(container, meta_entries_new, this->block, this->permutator.permutation_array);

	// Add meta records to the block buffers
	for(U32 i = 0; i < container.sizeWithoutCarryOver(); ++i) this->block += meta_entries_new[i];

	for(std::size_t i = 0; i < container.sizeWithoutCarryOver(); ++i) (meta_entries_new + i)->~MetaEntry();
	::operator delete[](static_cast<void*>(meta_entries_new));

	return true;
}

bool VariantImporter::AddRecord(const vcf_container_type& container, const U32 position, meta_type& meta){
	if(container.at(position)->pos > this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->n_bases){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->name << ':' << container.at(position)->pos + 1 << " > reported max size of contig (" << this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->n_bases + 1 << ")..." << std::endl;
		return false;
	}

	io::VcfGenotypeSummary s = container.GetGenotypeSummary(position, this->vcf_reader_->vcf_header_.GetNumberSamples());
	meta.controller.diploid          = (s.base_ploidy == 2);
	meta.controller.gt_mixed_phasing = s.mixed_phasing;
	meta.controller.gt_phase         = s.phase_if_uniform;
	meta.controller.mixed_ploidy     = (s.n_vector_end != 0);

	if(this->AddVcfFilterInfo(container.at(position), meta) == false) return false;
	if(this->AddVcfInfo(container.at(position), meta) == false) return false;
	if(this->AddVcfFormatInfo(container.at(position), meta) == false) return false;
	if(this->IndexRecord(container.at(position), meta) == false) return false;

	return true;
}

bool VariantImporter::AddVcfFilterInfo(const bcf1_t* record, meta_type& meta){
	// Add FILTER id list to the block. Filter information is unique in that the
	// data is not stored as (key,value)-tuples but as a key id. Because no data
	// is stored in the block, only the unique vectors of ids and their unique
	// ids are collected here. These keys are used to construct bit-vectors for
	// set-membership tests.
	std::vector<int> filter_ids;
	const int& n_filter_fields = record->d.n_flt;
	for(U32 i = 0; i < n_filter_fields; ++i){
		const int& hts_filter_key = record->d.flt[i]; // htslib IDX value
		const U32 global_key = this->filter_reorder_map_[hts_filter_key]; // tachyon global IDX value
		const U32 target_container = this->block.AddFilter(global_key);
		/*
		const reorder_map_type::const_iterator it = this->filter_local_map_.find(global_key); // search for global IDX
		if(it == this->filter_local_map_.end()){
			this->filter_local_map_[global_key] = this->filter_list_.size(); // local IDX
			this->filter_list_.push_back(global_key); // store local IDX at the key of the global IDX
			// Todo: spawn new object here if required.
		}
		const U32 target_container = this->filter_local_map_[global_key];
		*/
		assert(target_container < 65536);
		filter_ids.push_back(global_key);
	}

	return(this->AddVcfFilterPattern(filter_ids, meta));
}

bool VariantImporter::AddVcfInfo(const bcf1_t* record, meta_type& meta){
	// Add INFO fields to the block
	std::vector<int> info_ids;
	const int n_info_fields = record->n_info;
	for(U32 i = 0; i < n_info_fields; ++i){
		const int& hts_info_key = record->d.info[i].key; // htslib IDX value
		const U32 global_key = this->info_reorder_map_[hts_info_key]; // tachyon global IDX value
		const U32 target_container = this->block.AddInfo(global_key);
		/*
		const reorder_map_type::const_iterator it = this->info_local_map_.find(global_key); // search for global IDX
		if(it == this->info_local_map_.end()){
			this->info_local_map_[global_key] = this->info_list_.size(); // local IDX
			this->info_list_.push_back(global_key); // store local IDX at the key of the global IDX
		}
		const U32 target_container = this->info_local_map_[global_key];
		*/
		assert(target_container < 65536);
		info_ids.push_back(global_key);

		stream_container& destination_container = this->block.info_containers[target_container];
		const int& info_primitive_type = record->d.info[i].type;
		const int& stride_size         = record->d.info[i].len;
		const uint32_t& data_length    = record->d.info[i].vptr_len;
		const uint8_t* data            = record->d.info[i].vptr;
		int element_stride_size        = 0;

		//std::cerr << "Info:       " << hts_info_key << "(" << global_key << "," << (*this->block.info_map)[global_key] << "=" << target_container << ")->" << this->vcf_reader_->vcf_header_.GetInfo(hts_info_key)->id << " : " << io::BCF_TYPE_LOOKUP[record->d.info[i].type] << std::endl;
		//std::cerr << "Additional: " << stride_size << "," << data_length << " -> " << this->vcf_reader_->vcf_header_.GetNumberSamples() << std::endl;

		if(info_primitive_type == BCF_BT_INT8){
			element_stride_size = sizeof(int8_t);
			assert(data_length % element_stride_size == 0);
			const SBYTE* data_local = reinterpret_cast<const SBYTE*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
		} else if(info_primitive_type == BCF_BT_INT16){
			element_stride_size = sizeof(int16_t);
			assert(data_length % element_stride_size == 0);
			const S16* data_local = reinterpret_cast<const S16*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
		} else if(info_primitive_type == BCF_BT_INT32){
			element_stride_size = sizeof(int32_t);
			assert(data_length % element_stride_size == 0);
			const S32* data_local = reinterpret_cast<const S32*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
		} else if(info_primitive_type == BCF_BT_FLOAT){
			element_stride_size = sizeof(float);
			assert(data_length % element_stride_size == 0);
			const float* data_local = reinterpret_cast<const float*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
		} else if(info_primitive_type == BCF_BT_CHAR){
			element_stride_size = sizeof(char);
			const char* data_local = reinterpret_cast<const char*>(data);
			destination_container.AddCharacter(data_local, data_length);
		} else if(info_primitive_type == BCF_BT_NULL){
			element_stride_size = 0;
			assert(data_length == 0 && stride_size == 0);
		} else {
			std::cerr << utility::timestamp("ERROR","VCF") << "Unknown case: " << (int)info_primitive_type << std::endl;
			exit(1);
		}

		++destination_container;
		destination_container.AddStride(stride_size);
	}

	return(this->AddVcfInfoPattern(info_ids, meta));
}

bool VariantImporter::AddVcfFormatInfo(const bcf1_t* record, meta_type& meta){
	// Add FORMAT fields to the block
	std::vector<int> format_ids;
	const int n_format_fields = record->n_fmt;
	for(U32 i = 0; i < n_format_fields; ++i){
		const int& hts_format_key = record->d.fmt[i].id;; // htslib IDX value
		const U32 global_key = this->format_reorder_map_[hts_format_key]; // tachyon global IDX value
		const U32 target_container = this->block.AddFormat(global_key);
		/*
		const reorder_map_type::const_iterator it = this->format_local_map_.find(global_key); // search for global IDX
		if(it == this->format_local_map_.end()){
			this->format_local_map_[global_key] = this->format_list_.size(); // local IDX
			this->format_list_.push_back(global_key); // store local IDX at the key of the global IDX
		}
		const U32 target_container = this->format_local_map_[global_key];
		*/
		assert(target_container < 65536);
		format_ids.push_back(global_key);

		// Genotypes are a special case and are treated completely differently.
		// Because of this we simply skip that data here if it is available.
		if(this->vcf_reader_->vcf_header_.GetFormat(hts_format_key)->id == "GT"){
			continue;
		}

		stream_container& destination_container = this->block.format_containers[target_container];

		const int& format_primitive_type = record->d.fmt[i].type;
		const int& stride_size           = record->d.fmt[i].n;
		const uint32_t& data_length      = record->d.fmt[i].p_len;
		const uint8_t* data              = record->d.fmt[i].p;
		int element_stride_size          = 0;

		//std::cerr << "Format: " << hts_format_key << "(" << global_key << "," << this->format_local_map_[global_key] << "=" << target_container << ")->" << this->vcf_reader_->vcf_header_.GetFormat(hts_format_key)->id << " : " << io::BCF_TYPE_LOOKUP[record->d.fmt[i].type]
		//<< " stride size: " << stride_size << " data length: " << data_length << "->" << data_length/stride_size << std::endl;																																																			  //stream_container& target_container = this->block.format_containers[map_id];

		if(format_primitive_type == BCF_BT_INT8){
			element_stride_size = sizeof(int8_t);
			assert(data_length % element_stride_size == 0);
			const SBYTE* data_local = reinterpret_cast<const SBYTE*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
			assert(stride_size * element_stride_size * this->vcf_reader_->vcf_header_.GetNumberSamples() == data_length);
		} else if(format_primitive_type == BCF_BT_INT16){
			element_stride_size = sizeof(int16_t);
			assert(data_length % element_stride_size == 0);
			const S16* data_local = reinterpret_cast<const S16*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
			assert(stride_size * element_stride_size * this->vcf_reader_->vcf_header_.GetNumberSamples() == data_length);
		} else if(format_primitive_type == BCF_BT_INT32){
			element_stride_size = sizeof(int32_t);
			assert(data_length % element_stride_size == 0);
			const S32* data_local = reinterpret_cast<const S32*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
			assert(stride_size * element_stride_size * this->vcf_reader_->vcf_header_.GetNumberSamples() == data_length);
		} else if(format_primitive_type == BCF_BT_FLOAT){
			element_stride_size = sizeof(float);
			assert(data_length % element_stride_size == 0);
			const float* data_local = reinterpret_cast<const float*>(data);
			for(U32 j = 0; j < data_length/element_stride_size; ++j)
				destination_container.Add(data_local[j]);
			assert(stride_size * element_stride_size * this->vcf_reader_->vcf_header_.GetNumberSamples() == data_length);
		} else if(format_primitive_type == BCF_BT_CHAR){
			element_stride_size = sizeof(char);
			const char* data_local = reinterpret_cast<const char*>(data);
			destination_container.AddCharacter(data_local, data_length);
		} else {
			std::cerr << utility::timestamp("ERROR","VCF") << "Unknown case: " << (int)format_primitive_type << std::endl;
			exit(1);
		}

		++destination_container;
		destination_container.AddStride(stride_size);
	}

	return(this->AddVcfFormatPattern(format_ids, meta));
}

bool VariantImporter::AddVcfInfoPattern(const std::vector<int>& pattern, meta_type& meta){
	if(pattern.size()){
		meta.info_pattern_id = this->block.AddInfoPattern(pattern);
	}
	return true;
}

bool VariantImporter::AddVcfFormatPattern(const std::vector<int>& pattern, meta_type& meta){
	if(pattern.size()){
		meta.format_pattern_id = this->block.AddFormatPattern(pattern);
	}
	return true;
}

bool VariantImporter::AddVcfFilterPattern(const std::vector<int>& pattern, meta_type& meta){
	if(pattern.size()){
		meta.filter_pattern_id = this->block.AddFilterPattern(pattern);
	}
	return true;
}

bool VariantImporter::IndexRecord(const bcf1_t* record, const meta_type& meta){
	S32 index_bin = -1;

	// Ascertain that the meta entry has been evaluated
	// prior to executing this function.
	if(meta.n_alleles == 0){
		std::cerr << utility::timestamp("ERROR","IMPORT") << "The target meta record must be parsed prior to executing indexing functions..." << std::endl;
		return false;
	}

	S64 end_position_used = record->pos;

	// The Info field END is used as the end position of an internal if it is available. This field
	// is usually only set for non-standard variants such as SVs or other special meaning records.
	if(this->settings_.info_end_key != -1){
		// Linear search for the END key: this is not optimal but is probably faster
		// than first constructing a hash table for each record.
		const int n_info_fields = record->n_info;
		for(U32 i = 0; i < n_info_fields; ++i){
			if(record->d.info[i].key == this->settings_.info_end_key){
				U32 end = 0;
				switch(record->d.info[i].type){
				case(BCF_BT_INT8):  end = *reinterpret_cast<int8_t*> (record->d.info[i].vptr); break;
				case(BCF_BT_INT16): end = *reinterpret_cast<int16_t*>(record->d.info[i].vptr); break;
				case(BCF_BT_INT32): end = *reinterpret_cast<int32_t*>(record->d.info[i].vptr); break;
				default:
					std::cerr << "Illegal END primitive type: " << io::BCF_TYPE_LOOKUP[record->d.info[i].type] << std::endl;
					return false;
				}
				//std::cerr << "Found END at " << i << ".  END=" << end << " POS=" << record->pos + 1 << std::endl;
				index_bin = this->writer->index.index_[meta.contigID].add(record->pos, end, (U32)this->writer->index.current_block_number());
				//index_bin = 0;
				break;
			}
		}
	}

	// If the END field cannot be found then we check if the variant is a
	if(index_bin == -1){
		S32 longest = -1;
		// Iterate over available allele information and find the longest
		// SNV/indel length. The regex pattern ^[ATGC]{1,}$ searches for
		// simple SNV/indels.
		for(U32 i = 0; i < meta.n_alleles; ++i){
			if(std::regex_match(meta.alleles[i].allele, utility::YON_VARIANT_STANDARD)){
				if(meta.alleles[i].l_allele > longest)
					longest = meta.alleles[i].l_allele;
			}
		}

		// Update the variant index with the target bin(s) found.
		if(longest > 1){
			index_bin = this->writer->index.index_[meta.contigID].add(record->pos, record->pos + longest, (U32)this->writer->index.current_block_number());
			//index_bin = 0;
			end_position_used = record->pos + longest;
		}
		// In the cases of special-meaning alleles such as copy-number (e.g. <CN>)
		// or SV (e.g. A[B)) they are index according to their left-most value only.
		// This has the implication that they cannot be found by means of interval
		// intersection searches. If special-meaning variants were to be supproted
		// in the index then many more blocks would have to be searched for each
		// query as the few will dominate the many.
		else {
			index_bin = this->writer->index.index_[meta.contigID].add(record->pos, record->pos, (U32)this->writer->index.current_block_number());
			//index_bin = 0;
			//std::cerr << "fallback: " << record->pos+1 << std::endl;
		}

		//std::cerr << "End pos used: " << end_position_used << std::endl;
	}
	if(index_bin > this->index_entry.maxBin) this->index_entry.maxBin = index_bin;
	if(index_bin < this->index_entry.minBin) this->index_entry.minBin = index_bin;
	if(end_position_used > this->index_entry.maxPosition)
		this->index_entry.maxPosition = end_position_used;

	// Update number of entries in block
	++this->index_entry.n_variants;

	return true;
}

bool VariantImporter::addSite(meta_type& meta, bcf_entry_type& entry){
	// Assert position is in range
	if(entry.body->POS + 1 > this->header->getContig(entry.body->CHROM).bp_length){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << this->header->getContig(entry.body->CHROM).name << ':' << entry.body->POS+1 << " > reported max size of contig (" << this->header->getContig(entry.body->CHROM).bp_length << ")..." << std::endl;
		return false;
	}

	if(!this->parseBCFBody(meta, entry)){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Failed to encode BCF body..." << std::endl;
		return false;
	}

	// Add meta
	//this->block += meta;

	// Has END position?
	S32 index_bin = -1;

	S64 end_position_used = entry.body->POS;
	if(this->settings_.info_end_key != -1){
		// Linear search: this is not optimal but probably still faster
		// than generating a new hash table for each record
		for(U32 i = 0; i < entry.infoPointer; ++i){
			if(entry.infoID[i].mapID == this->settings_.info_end_key){
				const U32 end = entry.getInteger(entry.infoID[i].primitive_type, entry.infoID[i].l_offset);
				end_position_used = end;
				//std::cerr << "Found END at " << i << ".  END=" << end << " POS=" << entry.body->POS+1 << " BIN=" << reg2bin(entry.body->POS, end)  << std::endl;
				index_bin = this->writer->index.index_[meta.contigID].add(entry.body->POS, end, (U32)this->writer->index.current_block_number());
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
				if(meta.alleles[i].l_allele > longest) longest = meta.alleles[i].l_allele;
			} else {
				//std::cerr << "POS=" << entry.body->POS+1 << " no regex match: " << meta.alleles[i].toString() << std::endl;
			}
		}

		if(longest > 1){
			index_bin = this->writer->index.index_[meta.contigID].add(entry.body->POS, entry.body->POS + longest, (U32)this->writer->index.current_block_number());
			end_position_used = entry.body->POS + longest;
		} else { // fallback if all others fail
			index_bin = this->writer->index.index_[meta.contigID].add(entry.body->POS, entry.body->POS, (U32)this->writer->index.current_block_number());
		}
	}
	if(index_bin > this->index_entry.maxBin) this->index_entry.maxBin = index_bin;
	if(index_bin < this->index_entry.minBin) this->index_entry.minBin = index_bin;
	if(end_position_used > this->index_entry.maxPosition)
		this->index_entry.maxPosition = end_position_used;

	// Update number of entries in block
	++this->index_entry.n_variants;

	return true;
}

bool VariantImporter::parseBCFBody(meta_type& meta, bcf_entry_type& entry){
	for(U32 i = 0; i < entry.filterPointer; ++i){
		assert(entry.filterID[i].mapID != -1);
		this->block.addFieldFILTER(this->header->filter_remap[entry.filterID[i].mapID]);
	}

	for(U32 i = 0; i < entry.infoPointer; ++i){
		assert(entry.infoID[i].mapID != -1);
		const U32 mapID = this->block.addFieldINFO(this->header->info_remap[entry.infoID[i].mapID]);

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
		target_container.AddStride(entry.infoID[i].l_stride);
	}

	for(U32 i = 0; i < entry.formatPointer; ++i){
		assert(entry.formatID[i].mapID != -1);

		//const U32 mapID = this->block.format_fields.setGet(this->header->format_remap[entry.formatID[i].mapID]);
		const U32 mapID = this->block.addFieldFORMAT(this->header->format_remap[entry.formatID[i].mapID]);
		U32 internal_pos = entry.formatID[i].l_offset;

		// First value is always genotypes if there are any
		if(entry.hasGenotypes == true && i == 0)
			continue;

		// Hash INFO values
		stream_container& target_container = this->block.format_containers[mapID];

		// Flags and integers
		// These are BCF value types
		if(entry.formatID[i].primitive_type <= 3){
			for(U32 s = 0; s < this->header->n_samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container.Add(entry.getInteger(entry.formatID[i].primitive_type, internal_pos));
			}
		}
		// Floats
		else if(entry.formatID[i].primitive_type == bcf::BCF_FLOAT){
			for(U32 s = 0; s < this->header->n_samples; ++s){
				for(U32 j = 0; j < entry.formatID[i].l_stride; ++j)
					target_container.Add(entry.getFloat(internal_pos));
			}
		}
		// Chars
		else if(entry.formatID[i].primitive_type == bcf::BCF_CHAR){
			for(U32 s = 0; s < this->header->n_samples; ++s){
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
		target_container.AddStride(entry.formatID[i].l_stride);
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
