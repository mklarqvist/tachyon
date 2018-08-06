#include <fstream>
#include <regex>

#include "variant_importer.h"
#include "containers/checksum_container.h"

namespace tachyon {

VariantImporter::VariantImporter(const settings_type& settings) :
		settings_(settings),
		GT_available_(false),
		writer(nullptr)
{

}

VariantImporter::~VariantImporter(){
	delete this->writer;
}

void VariantImporter::clear(){
	this->vcf_container_.clear();
}

bool VariantImporter::Build(){
	if(!this->BuildVCF()){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
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

	for(U32 i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i){
		if(this->vcf_reader_->vcf_header_.contigs_[i].n_bases == 0){
			std::cerr << "No length declared for contig. Add maximum." << std::endl;
			this->vcf_reader_->vcf_header_.contigs_[i].n_bases = std::numeric_limits<S32>::max();
		}
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

	// Open a file handle or standard out for writing.
	if(!this->writer->open(this->settings_.output_prefix)){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	// Setup the encryption container.
	encryption::EncryptionDecorator encryption_manager;
	encryption::Keychain<> keychain;

	// Setup the checksums container.
	algorithm::VariantDigestManager checksums(YON_BLK_N_STATIC   + 1, // Add one for global checksum.
			this->vcf_reader_->vcf_header_.info_fields_.size()   + 1,
			this->vcf_reader_->vcf_header_.format_fields_.size() + 1);

	// The index needs to know how many contigs that's described in the
	// Vcf header and their lenghts. This information is needed to construct
	// the linear and quad-tree index most appropriate for the data.
	this->writer->index.Add(this->vcf_reader_->vcf_header_.contigs_);

	// Write out a fresh Tachyon header with the data from the Vcf header. As
	// this data will not be modified during the import stage it is safe to
	// write out now.
	this->writer->stream->write(&constants::FILE_HEADER[0], constants::FILE_HEADER_LENGTH); // Todo: fix
	this->WriteYonHeader();

	// Setup genotype permuter and genotype encoder.
	this->permutator.SetSamples(this->vcf_reader_->vcf_header_.GetNumberSamples());
	this->encoder.SetSamples(this->vcf_reader_->vcf_header_.GetNumberSamples());

	// Allocate containers and offsets for this file.
	// This is not strictly necessary but prevents nasty resize
	// calls in most cases.
	this->block.Allocate(this->vcf_reader_->vcf_header_.info_fields_.size(),
	                     this->vcf_reader_->vcf_header_.format_fields_.size(),
	                     this->vcf_reader_->vcf_header_.filter_fields_.size());

	// Resize containers
	const U32 resize_to = this->settings_.checkpoint_n_snps * sizeof(U32) * 2; // small initial allocation
	this->block.resize(resize_to);

	// Start porgress timer.
	algorithm::Timer timer; timer.Start();

	// Iterate over all available variants in the file or until encountering
	// an error.
	while(true){
		// Retrieve bcf1_t records using htslib and lazy evaluate them. Stop
		// after retrieving a set number of variants or if the interval between
		// the smallest and largest variant exceeds some distance in base pairs.
		if(this->vcf_container_.getVariants(this->settings_.checkpoint_n_snps,
		                                    this->settings_.checkpoint_bases,
		                                    this->vcf_reader_) == false){
			break;
		}

		// This pointer here is borrowed from the PPA manager
		// during import stages. Do not destroy the target block
		// before finishing with this.
		this->block.gt_ppa = &this->permutator.permutation_array;

		if(this->GT_available_ && this->settings_.permute_genotypes){
			// Only store the permutation array if the number of samples
			// are greater then one (1).
			if(this->vcf_reader_->vcf_header_.GetNumberSamples() > 1){
				if(this->permutator.Build(this->vcf_container_, this->vcf_reader_->vcf_header_) == false)
					return false;

				this->block.header.controller.hasGTPermuted = true;
			}
		}

		if(this->AddRecords(this->vcf_container_) == false) return false;

		this->block.header.controller.hasGT  = this->GT_available_;
		this->block.header.n_variants        = this->vcf_container_.sizeWithoutCarryOver();
		this->block.UpdateContainers();
		this->block.Finalize();

		// Perform compression using standard parameters.
		if(!this->compression_manager.Compress(this->block, this->settings_.compression_level, 6)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Encrypt the variant block if desired.
		if(this->settings_.encrypt_data){
			// Generate field-unique identifiers.
			this->GenerateIdentifiers();

			// Start encryption.
			this->block.header.controller.anyEncrypted = true;
			if(!encryption_manager.encrypt(this->block, keychain, YON_ENCRYPTION_AES_256_GCM)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to encrypt..." << std::endl;
			}
		}

		// Write the current variant block.
		this->WriteBlock();

		// Update checksums container with the available data.
		checksums += this->block;

		if(!SILENT){
			std::cerr << utility::timestamp("PROGRESS") <<
			std::setfill(' ') << std::setw(10) << this->writer->n_variants_written << ' ' <<
			std::setfill(' ') << std::setw(10) << utility::toPrettyDiskString(this->writer->stream->tellp()) << '\t' <<
			timer.ElapsedString() << ' ' <<
			this->vcf_reader_->vcf_header_.GetContig(this->vcf_container_.front()->rid)->name << ":" << this->vcf_container_.front()->pos + 1 << "->" << this->vcf_container_.back()->pos + 1 << std::endl;
		}

		// Clear current data.
		this->clear();
		this->block.clear();
		this->index_entry.reset();
	}
	// Do not delete the borrowed pointer.
	this->block.gt_ppa = nullptr;

	// Finalize writing procedure.
	this->WriteFinal(checksums);
	this->WriteKeychain(keychain);

	// All done
	return(true);
}

bool VariantImporter::AddRecords(const vcf_container_type& container){
	meta_type* meta_entries = static_cast<meta_type*>(::operator new[](container.sizeWithoutCarryOver() * sizeof(meta_type)));
	for(U32 i = 0; i < container.sizeWithoutCarryOver(); ++i){
		// Transmute a bcf record into a meta structure
		new( meta_entries + i ) meta_type( this->vcf_container_[i], this->block.header.minPosition );

		// Add the record data
		if(this->AddRecord(container, i, meta_entries[i]) == false)
			return false;
	}

	// Add FORMAT:GT field data.
	this->encoder.Encode(container, meta_entries, this->block, this->permutator.permutation_array);

	// Add meta records to the block buffers
	for(U32 i = 0; i < container.sizeWithoutCarryOver(); ++i) this->block += meta_entries[i];

	for(std::size_t i = 0; i < container.sizeWithoutCarryOver(); ++i) (meta_entries + i)->~MetaEntry();
	::operator delete[](static_cast<void*>(meta_entries));

	return true;
}

bool VariantImporter::AddRecord(const vcf_container_type& container, const U32 position, meta_type& meta){
	if(container.at(position)->pos > this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->n_bases){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->name << ':' << container.at(position)->pos + 1 << " > reported max size of contig (" << this->vcf_reader_->vcf_header_.GetContig(container.at(position)->rid)->n_bases + 1 << ")..." << std::endl;
		return false;
	}

	if(this->AddVcfFilterInfo(container.at(position), meta) == false) return false;
	if(this->AddVcfInfo(container.at(position), meta) == false) return false;

	if(container.at(position)->n_fmt){
		if(this->AddVcfFormatInfo(container.at(position), meta) == false) return false;

		// Perform these actions if FORMAT:GT data is available.
		const int& hts_format_key = container.at(position)->d.fmt[0].id; // htslib IDX value
		if(this->vcf_reader_->vcf_header_.GetFormat(hts_format_key)->id != "GT"){
			meta.controller.gt_available = false;
		} else
			meta.controller.gt_available = true;
	}

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
		assert(target_container < 65536);
		info_ids.push_back(global_key);

		stream_container& destination_container = this->block.info_containers[target_container];
		const int& info_primitive_type = record->d.info[i].type;
		const int& stride_size         = record->d.info[i].len;
		const uint32_t& data_length    = record->d.info[i].vptr_len;
		const uint8_t* data            = record->d.info[i].vptr;
		int element_stride_size        = 0;

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
			index_bin = this->writer->index.index_[meta.contigID].add(record->pos,
			                                                          record->pos + longest,
			                                                          (U32)this->writer->index.current_block_number());
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

bool VariantImporter::UpdateIndex(){
	this->index_entry.blockID         = this->block.header.block_hash;
	this->index_entry.byte_offset_end = this->writer->stream->tellp();
	this->index_entry.contigID        = this->vcf_container_.front()->rid;
	this->index_entry.minPosition     = this->vcf_container_.front()->pos;
	this->index_entry.n_variants      = this->vcf_container_.sizeWithoutCarryOver();
	this->writer->index.index_.linear_at(index_entry.contigID) += this->index_entry;
	this->index_entry.reset();
	++this->writer->n_blocks_written;
	this->writer->n_variants_written += this->vcf_container_.sizeWithoutCarryOver();
	++this->writer->index.number_blocks;

	return true;
}

bool VariantImporter::WriteBlock(){
	this->index_entry.byte_offset = this->writer->stream->tellp();
	this->block.write(*this->writer->stream, this->stats_basic, this->stats_info, this->stats_format);

	// After all compression and writing is finished the header
	// offsets are themselves compressed and stored in the block.
	this->block.PackFooter(); // Pack footer into buffer.
	this->compression_manager.zstd_codec.Compress(block.footer_support);
	this->writer->WriteBlockFooter(this->block.footer_support);
	this->writer->WriteEndOfBlock(); // End-of-block marker
	this->UpdateIndex(); // Update index.
	return(this->writer->stream->good());
}

bool VariantImporter::WriteFinal(algorithm::VariantDigestManager& checksums){
	// Done importing
	this->writer->stream->flush();

	core::Footer footer;
	footer.offset_end_of_data = this->writer->stream->tellp();
	footer.n_blocks           = this->writer->n_blocks_written;
	footer.n_variants         = this->writer->n_variants_written;
	assert(footer.n_blocks == this->writer->index.GetLinearSize());

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
	return(this->writer->stream->good());
}

bool VariantImporter::WriteKeychain(const encryption::Keychain<>& keychain){
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
	return true;
}

bool VariantImporter::WriteYonHeader(){
	VariantHeader yon_header(this->vcf_reader_->vcf_header_);

	io::BasicBuffer temp(500000);
	io::BasicBuffer temp_cmp(temp);
	temp << yon_header;
	this->compression_manager.zstd_codec.Compress(temp, temp_cmp, 20);
	uint32_t l_data   = temp.size();
	uint32_t l_c_data = temp_cmp.size();
	utility::SerializePrimitive(l_data,   *this->writer->stream);
	utility::SerializePrimitive(l_c_data, *this->writer->stream);
	this->writer->stream->write(temp_cmp.data(), l_c_data);
	return(this->writer->stream->good());
}

bool VariantImporter::GenerateIdentifiers(void){
	BYTE RANDOM_BYTES[32];
	for(U32 i = 0; i < this->vcf_container_.sizeWithoutCarryOver(); ++i){
		U64 b_hash;
		while(true){
			RAND_bytes(&RANDOM_BYTES[0], 32);
			b_hash = XXH64(&RANDOM_BYTES[0], 32, 1337);
			hash_map_type::const_iterator it = this->block_hash_map.find(b_hash);
			if(it == this->block_hash_map.end()){
				this->block_hash_map[b_hash] = 0; // Number doesn't matter.
				break;
			}
		}
	}
	return true;
}

} /* namespace Tachyon */
