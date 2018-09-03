#include "variant_reader.h"


namespace tachyon{

VariantReader::VariantReader() :
	b_data_start(0)
{}

VariantReader::VariantReader(const std::string& filename) :
	b_data_start(0),
	basic_reader(filename)
{}

VariantReader::~VariantReader(){}

VariantReader::VariantReader(const self_type& other) :
	b_data_start(other.b_data_start),
	basic_reader(other.basic_reader),
	//variant_container(other.variant_container),

	block_settings(other.block_settings),
	settings(other.settings),
	//variant_filters(other.variant_filters), // illeal to copy

	global_header(other.global_header),
	global_footer(other.global_footer),
	index(other.index),
	checksums(other.checksums),
	keychain(other.keychain),
	interval_container(other.interval_container),
	occ_table(other.occ_table)
{
	this->basic_reader.stream_.seekg(this->b_data_start);
	this->variant_container << this->global_header;
}

bool VariantReader::open(void){
	if(this->settings.input.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	if(this->basic_reader.open() == false){
		std::cerr << "Failed to open" << std::endl;
		return false;
	}

	if(this->basic_reader.filesize_ <= YON_FOOTER_LENGTH){
		std::cerr << utility::timestamp("ERROR") << "File is corrupted!" << std::endl;
		return false;
	}

	// Seek to start of footer
	this->basic_reader.stream_.seekg((uint64_t)this->basic_reader.filesize_ - YON_FOOTER_LENGTH);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to seek in file!" << std::endl;
		return false;
	}
	this->basic_reader.stream_ >> this->global_footer;

	// Validate footer
	if(this->global_footer.validate() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate footer!" << std::endl;
		return false;
	}

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to read file!" << std::endl;
		return false;
	}

	// Seek to start of file
	this->basic_reader.stream_.seekg(0);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to rewind file!" << std::endl;
		return false;
	}

	// Load header
	//this->stream >> this->global_header;
	char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	this->basic_reader.stream_.read(&magic_string[0], tachyon::constants::FILE_HEADER_LENGTH);
	if(strncmp(&magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate Tachyon magic string!" << std::endl;
		return false;
	}

	uint32_t l_data   = 0;
	uint32_t l_c_data = 0;
	utility::DeserializePrimitive(l_data, this->basic_reader.stream_);
	utility::DeserializePrimitive(l_c_data, this->basic_reader.stream_);

	io::BasicBuffer header_uncompressed(l_data + 1024);
	io::BasicBuffer header_compressed(l_c_data + 1024); header_compressed.n_chars_ = l_c_data;

	this->basic_reader.stream_.read(header_compressed.data(), l_c_data);

	if(!this->codec_manager.zstd_codec.Decompress(header_compressed, header_uncompressed)){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	assert(header_uncompressed.size() == l_data);
	header_uncompressed >> this->global_header; // parse header from buffer

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to get header!" << std::endl;
		return false;
	}

	// Keep track of start of data offset in the byte stream.
	// We use this information if spawning copies of this
	// reader object.
	this->b_data_start = this->basic_reader.stream_.tellg();

	// Add global header pointer to variant container.
	this->variant_container << this->global_header;

	// Keep track of start position
	const uint64_t return_pos = this->basic_reader.stream_.tellg();
	this->basic_reader.stream_.seekg(this->global_footer.offset_end_of_data);
	this->basic_reader.stream_ >> this->index;
	this->basic_reader.stream_ >> this->checksums;
	this->basic_reader.stream_.seekg(return_pos);

	return(this->basic_reader.stream_.good());
}

bool VariantReader::NextBlock(){
	if(this->CheckNextValid() == false) return false;

	// Reset and re-use
	this->variant_container.reset();

	if(!this->variant_container.GetBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return false;

	if(!this->codec_manager.zstd_codec.Decompress(this->variant_container.GetBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
		return false;
	}
	this->variant_container.GetBlock().footer_support.data_uncompressed >> this->variant_container.GetBlock().footer;

	// Attempts to read a YON block with the settings provided
	if(!this->variant_container.ReadBlock(this->basic_reader.stream_, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.AnyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.Decrypt(this->variant_container.GetBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.Decompress(this->variant_container.GetBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

containers::VariantBlockContainer VariantReader::ReturnBlock(void){
	variant_container_type c;
	c << this->GetGlobalHeader();

	if(this->CheckNextValid() == false)
		return c;

	if(!c.GetBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return(c);

	if(!this->codec_manager.zstd_codec.Decompress(c.GetBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	c.GetBlock().footer_support.data_uncompressed >> c.GetBlock().footer;

	// Attempts to read a YON block with the settings provided
	if(!c.ReadBlock(this->basic_reader.stream_, this->block_settings))
		return(c);

	// encryption manager ascertainment
	if(c.AnyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return(c);
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.Decrypt(c.GetBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return(c);
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.Decompress(c.GetBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return(c);
	}

	return(c);
}

bool VariantReader::GetBlock(const index_entry_type& index_entry){
	// If the stream is not good then return.
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// Seek to target block id with the help of the linear index.
	this->basic_reader.stream_.seekg(index_entry.byte_offset);
	if(this->basic_reader.stream_.good() == false){
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to seek to given offset using target index entry!" << std::endl;
		return(false);
	}

	// Load the next block if possible.
	return(this->NextBlock());
}

bool VariantReader::LoadKeychainFile(void){
	std::ifstream keychain_reader(settings.keychain_file, std::ios::binary | std::ios::in);
	if(!keychain_reader.good()){
		std::cerr << tachyon::utility::timestamp("ERROR") <<  "Failed to open keychain: " << settings.keychain_file << "..." << std::endl;
		return false;
	}

	keychain_reader >> this->keychain;
	if(!keychain_reader.good()){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse keychain..." << std::endl;
		return false;
	}
	return true;
}

TACHYON_VARIANT_CLASSIFICATION_TYPE VariantReader::ClassifyVariant(const meta_entry_type& meta, const uint32_t& allele) const{
	const int32_t ref_size = meta.alleles[0].size();
	const int32_t l_diff   = ref_size - meta.alleles[allele].size();

	if(meta.alleles[0].allele[0] == '<' || meta.alleles[allele].allele[0] == '<')
		return(YON_VARIANT_CLASS_SV);
	else if(l_diff == 0){
		if(ref_size == 1 && meta.alleles[0].allele[0] != meta.alleles[allele].allele[0]){
			if(meta.alleles[allele].allele[0] == 'A' ||
			   meta.alleles[allele].allele[0] == 'T' ||
			   meta.alleles[allele].allele[0] == 'G' ||
			   meta.alleles[allele].allele[0] == 'C')
			{
				return(YON_VARIANT_CLASS_SNP);
			}
			else return(YON_VARIANT_CLASS_UNKNOWN);
		}
		else if(ref_size != 1){
			uint32_t n_characters_identical = 0;
			const uint32_t length_shortest = ref_size < meta.alleles[allele].size()
					                    ? ref_size
					                    : meta.alleles[allele].size();

			for(uint32_t c = 0; c < length_shortest; ++c)
				n_characters_identical += (meta.alleles[0].allele[c] == meta.alleles[allele].allele[c]);

			if(n_characters_identical == 0) return(YON_VARIANT_CLASS_MNP);
			else return(YON_VARIANT_CLASS_CLUMPED);
		}
	} else {
		const uint32_t length_shortest = ref_size < meta.alleles[allele].size()
		                            ? ref_size
		                            : meta.alleles[allele].size();

		// Keep track of non-standard characters.
		uint32_t n_characters_non_standard = 0;

		// Iterate over available characters and check for non-standard
		// genetic characters (ATGC).
		for(uint32_t c = 0; c < length_shortest; ++c){
			n_characters_non_standard += (meta.alleles[allele].allele[c] != 'A' &&
			                              meta.alleles[allele].allele[c] != 'T' &&
			                              meta.alleles[allele].allele[c] != 'C' &&
			                              meta.alleles[allele].allele[c] != 'G');
		}

		// If non-standard characters are found then return as unknown
		// type. Otherwise, return classification as an indel.
		if(n_characters_non_standard) return(YON_VARIANT_CLASS_UNKNOWN);
		else return(YON_VARIANT_CLASS_INDEL);
	}
	return(YON_VARIANT_CLASS_UNKNOWN);
}

void VariantReader::OuputVcfWrapper(io::BasicBuffer& output_buffer, yon1_t& entry) const{
	utility::ToVcfString(output_buffer, '\t', *entry.meta, this->global_header);
	output_buffer += '\t';

	// Print Filter, Info, and Format if available.
	this->OutputFilterVcf(output_buffer, entry);
	this->OutputInfoVcf(output_buffer, entry);
	this->OutputFormatVcf(output_buffer, entry);
	output_buffer += '\n';

	if(output_buffer.size() > 65536){
		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
	}
}

void VariantReader::OutputInfoVcf(io::BasicBuffer& output_buffer, yon1_t& entry) const{
	// Print Info.
	if(entry.n_info){
		const uint32_t n_info_avail = entry.info_ids->size();
		if(n_info_avail){
			if(entry.info_hdr[0]->yon_type == YON_VCF_HEADER_FLAG){
				output_buffer += entry.info_hdr[0]->id;
			} else {
				output_buffer += entry.info_hdr[0]->id;
				output_buffer += '=';
				entry.info[0]->ToVcfString(output_buffer);
			}

			for(uint32_t j = 1; j < n_info_avail; ++j){
				output_buffer += ';';
				if(entry.info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
					output_buffer += entry.info_hdr[j]->id;
				} else {
					output_buffer += entry.info_hdr[j]->id;
					output_buffer += '=';
					entry.info[j]->ToVcfString(output_buffer);
				}
			}

			if(this->GetBlockSettings().annotate_extra){
				output_buffer += ';';
				entry.EvaluateSummary(true);
				entry.gt_sum->d->PrintVcf(output_buffer);
			}
		}
	} else {
		if(this->GetBlockSettings().annotate_extra){
			entry.EvaluateSummary(true);
			entry.gt_sum->d->PrintVcf(output_buffer);
		} else
			output_buffer += '.';
	}
}

void VariantReader::OutputFormatVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const{
	if(entry.n_format){
		output_buffer += '\t';

		const uint32_t n_format_avail = entry.format_ids->size();
		if(n_format_avail){
			output_buffer += entry.format_hdr[0]->id;
			for(uint32_t j = 1; j < n_format_avail; ++j){
				output_buffer += ':';
				output_buffer += entry.format_hdr[j]->id;
			}
			output_buffer += '\t';

			// Vcf FORMAT values are interleaved such that values are
			// presented in a columnar representation with data for
			// each sample is concatenated together such as the pattern
			// GT:AF:AS is display for three samples as:
			//
			// Sample 1    Sample 2    Sample 3
			// 0|0:0.5|ABC 0|0:0.5|ABC 0|0:0.5|ABC
			//
			// This memory layout requires the interleaving of the
			// internally separated data streams resulting in slower
			// Vcf printing speeds compared to the naive Bcf file format.
			//
			// First calculate the FORMAT:GT field for this variant site.
			// Case when the only available FORMAT field is the GT field.
			if(n_format_avail == 1 && entry.is_loaded_gt &&
			   entry.meta->controller.gt_available &&
			   (this->GetBlockSettings().display_static & YON_BLK_BV_GT))
			{
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				// Iterate over samples and print FORMAT:GT value in Vcf format.
				entry.gt->d_exp[0]->PrintVcf(output_buffer, entry.gt->m);
				for(uint32_t s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.gt->d_exp[s]->PrintVcf(output_buffer, entry.gt->m);
				}

				entry.gt->d_exp = nullptr;
			}
			// Case when there are > 1 Vcf Format fields and the GT field
			// is available.
			else if(n_format_avail > 1 && entry.is_loaded_gt &&
			        entry.meta->controller.gt_available &&
			        (this->GetBlockSettings().display_static & YON_BLK_BV_GT))
			{
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				entry.gt->d_exp[0]->PrintVcf(output_buffer, entry.gt->m);
				for(uint32_t g = 1; g < n_format_avail; ++g){
					output_buffer += ':';
					entry.fmt[g]->ToVcfString(output_buffer, 0);
				}
				for(uint32_t s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.gt->d_exp[s]->PrintVcf(output_buffer, entry.gt->m);
					for(uint32_t g = 1; g < n_format_avail; ++g){
						output_buffer += ':';
						entry.fmt[g]->ToVcfString(output_buffer, s);
					}
				}

				entry.gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				entry.fmt[0]->ToVcfString(output_buffer, 0);
				for(uint32_t g = 1; g < n_format_avail; ++g){
					output_buffer += ':';
					entry.fmt[g]->ToVcfString(output_buffer, 0);
				}

				for(uint32_t s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.fmt[0]->ToVcfString(output_buffer, s);
					for(uint32_t g = 1; g < n_format_avail; ++g){
						output_buffer += ':';
						entry.fmt[g]->ToVcfString(output_buffer, s);
					}
				}
			}
		}
	}
}

void VariantReader::OutputFilterVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const{
	if(entry.n_filter){
		const uint32_t n_filter_avail = entry.filter_ids->size();
		if(n_filter_avail){
			output_buffer += entry.filter_hdr[0]->id;
			for(uint32_t j = 1; j < n_filter_avail; ++j){
				output_buffer += ';';
				output_buffer += entry.filter_hdr[j]->id;
			}
		} else {
			output_buffer += '.';
		}
	} else output_buffer += '.';
	output_buffer += '\t';
}

uint64_t VariantReader::OutputVcfLinear(void){
	this->variant_container.AllocateGenotypeMemory();
	// temp
	//if(this->occ_table.ReadTable("/media/mdrk/NVMe/1kgp3/populations/integrated_call_samples_v3.20130502.ALL.panel", this->GetGlobalHeader(), '\t') == false){
	//	return(0);
	//}

	while(this->NextBlock()){
		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);
		io::BasicBuffer output_buffer(100000);
		// If occ table is built.
		//objects->occ = &occ;
		//objects->EvaluateOcc(this->GetCurrentContainer().GetBlock().gt_ppa);

		for(uint32_t i = 0; i < objects->meta_container->size(); ++i){
			if(this->variant_filters.Filter(entries[i], i) == false)
				continue;

			// Each entry evaluate occ if available.
			//entries[i].occ = objects->occ;
			//entries[i].EvaluateOcc();

			this->OuputVcfWrapper(output_buffer, entries[i]);
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		delete [] entries;
		//objects->occ = nullptr;
		delete objects;
	}

	return 0;
}

uint64_t VariantReader::OutputVcfSearch(void){
	this->variant_container.AllocateGenotypeMemory();

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;

	for(uint32_t i = 0; i < this->interval_container.GetBlockList().size(); ++i){
		this->GetBlock(this->interval_container.GetBlockList()[i]);

		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);
		io::BasicBuffer output_buffer(100000);

		for(uint32_t i = 0; i < objects->meta_container->size(); ++i){
			if((this->*filter_intervals)(objects->meta_container->at(i)) == false)
				continue;

			if(this->variant_filters.Filter(entries[i], i) == false)
				continue;

			this->OuputVcfWrapper(output_buffer, entries[i]);
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		delete [] entries;
		delete objects;
	}

	return 0;
}

uint64_t VariantReader::OutputRecords(void){
	this->UpdateHeaderView();
	this->interval_container.Build(this->global_header);

	// Load encryption keychain if available.
	if(this->settings.keychain_file.size())
		this->LoadKeychainFile();

	if(this->settings.use_htslib){
		if(this->interval_container.size()) return(this->OutputHtslibVcfSearch());
		else return(this->OutputHtslibVcfLinear());
	}

	if(this->GetBlockSettings().show_vcf_header)
		this->global_header.PrintVcfHeader(std::cout);

	if(this->interval_container.size()) return(this->OutputVcfSearch());
	else return(this->OutputVcfLinear());
}

uint64_t VariantReader::OutputHtslibVcfLinear(void){
	this->variant_container.AllocateGenotypeMemory();

	// Open a htslib file handle for the target output
	// destination.
	char hts_stream_type[2];
	hts_stream_type[0] = 'w'; hts_stream_type[1] = this->settings.output_type;
	htsFile *fp = hts_open(this->settings.output.c_str(), hts_stream_type);

	// Convert the internal yon header to a bcf_hdr_t
	// structure.
	bcf_hdr_t* hdr = this->GetGlobalHeader().ConvertVcfHeader(!this->settings.drop_format);
	if ( bcf_hdr_write(fp, hdr) != 0 ) {
		std::cerr << "Failed to write header to " << this->settings.output << std::endl;
		exit(1);
	}

	// Initialize an empty record that we will keep
	// reusing as we iterate over available yon records.
	bcf1_t *rec = bcf_init1();

	// Iterate over available blocks.
	while(this->NextBlock()){
		// Lazy evaluate yon records.
		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);

		// Iterate over available records in this block.
		for(uint32_t i = 0; i < objects->meta_container->size(); ++i){
			if(this->variant_filters.Filter(entries[i], i) == false)
				continue;

			entries[i].meta->UpdateHtslibVcfRecord(rec, hdr);
			this->OutputHtslibVcfInfo(rec, hdr, entries[i]);
			this->OutputHtslibVcfFormat(rec, hdr, entries[i]);
			this->OutputHtslibVcfFilter(rec, hdr, entries[i]);

			if ( bcf_write1(fp, hdr, rec) != 0 ){
				std::cerr << "Failed to write record to " << this->settings.output;
				exit(1);
			}

			bcf_clear1(rec);
		}

		// Cleanup lazy evaluation of yon records.
		delete [] entries;
		delete objects;
	}

	// Cleanup htslib bcf1_t and bcf_hdr_t structures.
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);

	// Close file handle.
	int ret;
	if ( (ret=hts_close(fp)) ) {
		fprintf(stderr,"hts_close(%s): non-zero status %d\n",this->settings.output.data(),ret);
		exit(ret);
	}

	return 0;
}

uint64_t VariantReader::OutputHtslibVcfSearch(void){
	this->variant_container.AllocateGenotypeMemory();

	// Open a htslib file handle for the target output
	// destination.
	char hts_stream_type[2];
	hts_stream_type[0] = 'w'; hts_stream_type[1] = this->settings.output_type;
	htsFile *fp = hts_open(this->settings.output.c_str(), hts_stream_type);

	// Convert the internal yon header to a bcf_hdr_t
	// structure.
	bcf_hdr_t* hdr = this->GetGlobalHeader().ConvertVcfHeader(!this->settings.drop_format);
	if ( bcf_hdr_write(fp, hdr) != 0 ) {
		std::cerr << "Failed to write header to " << this->settings.output << std::endl;
		exit(1);
	}

	// Initialize an empty record that we will keep
	// reusing as we iterate over available yon records.
	bcf1_t *rec = bcf_init1();

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;

	// Iterate over available blocks.
	while(this->NextBlock()){
		// Lazy evaluate yon records.
		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);

		// Iterate over available records in this block.
		for(uint32_t i = 0; i < objects->meta_container->size(); ++i){
			if((this->*filter_intervals)(objects->meta_container->at(i)) == false)
				continue;

			if(this->variant_filters.Filter(entries[i], i) == false)
				continue;

			entries[i].meta->UpdateHtslibVcfRecord(rec, hdr);
			this->OutputHtslibVcfInfo(rec, hdr, entries[i]);
			this->OutputHtslibVcfFormat(rec, hdr, entries[i]);
			this->OutputHtslibVcfFilter(rec, hdr, entries[i]);

			if ( bcf_write1(fp, hdr, rec) != 0 ){
				std::cerr << "Failed to write record to " << this->settings.output;
				exit(1);
			}

			bcf_clear1(rec);
		}

		// Cleanup lazy evaluation of yon records.
		delete [] entries;
		delete objects;
	}

	// Cleanup htslib bcf1_t and bcf_hdr_t structures.
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);

	// Close file handle.
	int ret;
	if ( (ret=hts_close(fp)) ) {
		fprintf(stderr,"hts_close(%s): non-zero status %d\n",this->settings.output.data(),ret);
		exit(ret);
	}

	return 0;
}

void VariantReader::OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr, yon1_t& entry) const{
	if(entry.n_info){
		const uint32_t n_info_avail = entry.info_ids->size();
		if(n_info_avail){
			for(uint32_t j = 0; j < n_info_avail; ++j){
				if(entry.info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
					bcf_update_info_flag(hdr, rec, entry.info_hdr[j]->id.data(), NULL, 1);
				} else {
					entry.info[j]->UpdateHtslibVcfRecordInfo(rec, hdr, entry.info_hdr[j]->id);
				}
			}

			if(this->GetBlockSettings().annotate_extra){
				entry.EvaluateSummary(true);
				entry.gt_sum->d->UpdateHtslibVcfRecord(rec, hdr);
			}
		}
	} else {
		if(this->GetBlockSettings().annotate_extra){
			entry.EvaluateSummary(true);
			entry.gt_sum->d->UpdateHtslibVcfRecord(rec, hdr);
		}
	}
}

void VariantReader::OutputHtslibVcfFormat(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_t& entry) const{
	if(entry.n_format){
		const uint32_t n_format_avail = entry.format_ids->size();
		if(n_format_avail){
			// Case when the only available FORMAT field is the GT field.
			if(n_format_avail == 1 && entry.is_loaded_gt &&
			   entry.meta->controller.gt_available &&
			   (this->GetBlockSettings().display_static & YON_BLK_BV_GT))
			{
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();
				entry.gt->UpdateHtslibGenotypes(rec, hdr);
				entry.gt->d_exp = nullptr;
			}
			// Case when there are > 1 Vcf Format fields and the GT field
			// is available.
			else if(n_format_avail > 1 && entry.is_loaded_gt &&
			        entry.meta->controller.gt_available &&
			        (this->GetBlockSettings().display_static & YON_BLK_BV_GT))
			{
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				entry.gt->UpdateHtslibGenotypes(rec, hdr);
				for(uint32_t g = 1; g < n_format_avail; ++g)
					entry.format_containers[g]->UpdateHtslibVcfRecord(entry.id_block, rec, hdr, entry.format_hdr[g]->id);

				entry.gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				for(uint32_t g = 0; g < n_format_avail; ++g)
					entry.format_containers[g]->UpdateHtslibVcfRecord(entry.id_block, rec, hdr, entry.format_hdr[g]->id);
			}
		}
	}
}

void VariantReader::OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_t& entry) const{
	if(entry.n_filter){
		for(uint32_t k = 0; k < entry.filter_ids->size(); ++k){
			int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, entry.filter_hdr[k]->id.data());
			bcf_update_filter(hdr, rec, &tmpi, 1);
		}
	}
}

void VariantReader::UpdateHeaderView(void){
	io::VcfExtra e;
	e.key = "tachyon_viewVersion";
	e.value = tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	this->GetGlobalHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetGlobalHeader().extra_fields_.push_back(e);
	e.key = "tachyon_viewCommand";
	e.value = tachyon::constants::LITERAL_COMMAND_LINE;
	this->GetGlobalHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetGlobalHeader().extra_fields_.push_back(e);
	e.key = "tachyon_viewCommandSettings";
	e.value = this->GetSettings().get_settings_string();
	this->GetGlobalHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetGlobalHeader().extra_fields_.push_back(e);
}

bool VariantReader::Stats(void){
	this->variant_container.AllocateGenotypeMemory();
	// temp
	//if(this->occ_table.ReadTable("/media/mdrk/NVMe/1kgp3/populations/integrated_call_samples_v3.20130502.ALL.panel", this->GetGlobalHeader(), '\t') == false){
	//	return(0);
	//}

	/*
	uint32_t n_blocks = 0;
	for(uint32_t i = 0; i < this->GetIndex().index_.n_contigs_; ++i){
		n_blocks += this->GetIndex().index_.linear_[i].size();
	}
	uint32_t n_threads = std::thread::hardware_concurrency();

	std::cerr << "blocks: " << n_blocks << " and threads " << n_threads << std::endl;

	std::vector<std::pair<uint32_t, uint32_t>> workload(n_threads);
	uint32_t n_cumulative = 0;
	assert(n_threads != 0);
	uint32_t n_step_size = n_blocks / n_threads;
	for(uint32_t i = 0; i < n_threads - 1; ++i){
		workload[i] = std::pair<uint32_t, uint32_t>(n_cumulative, n_cumulative + n_step_size);
		n_cumulative += n_step_size;
	}
	workload.back().first = n_cumulative;
	workload.back().second = n_blocks;

	for(uint32_t i = 0; i < n_threads; ++i){
		std::cerr << i << ": " << workload[i].first << "-" << workload[i].second << std::endl;
	}

	// Final solution.
	yon_stats_tstv s(this->GetGlobalHeader().GetNumberSamples());

	//
	yon_stats_thread_pool pool(n_threads);
	std::vector<std::thread*> threads(n_threads);
	for(uint32_t i = 0; i < n_threads; ++i){
		pool[i].reader      = new VariantReader(*this);
		pool[i].block_from  = workload[i].first;
		pool[i].block_to    = workload[i].second;
		threads[i] = pool[i].Start();
	}

	for(uint32_t i = 0; i < n_threads; ++i) threads[i]->join();
	std::cerr << "done joining" << std::endl;

	return true;
	*/
	yon_stats_tstv s(this->GetGlobalHeader().GetNumberSamples());

	// stupid test
	//VariantReader v(*this);

	while(this->NextBlock()){
		//variant_container_type c = v.ReturnBlock();
		//variant_container_type& c = this->GetCurrentContainer();
		//c.AllocateGenotypeMemory();

		// Move all data to new instance.
		//variant_container_type c(std::move(this->variant_container));
		//variant_container_type c;
		//c = std::move(this->variant_container);
		//variant_container_type c(this->variant_container);

		objects_type* objects = this->variant_container.LoadObjects(this->block_settings);
		yon1_t* entries = this->variant_container.LazyEvaluate(*objects);

		// Debug
		//std::cerr << objects->meta_container->front().position << "->" << objects->meta_container->back().position << std::endl;

		// Debug
		//containers::DataContainer dc2 = objects->format_containers[1]->ToDataContainer();
		//std::cerr << "fmt1\t" << dc2.header.n_additions << "," << dc2.header.n_strides << "," << dc2.GetSizeUncompressed() << std::endl;



		// If occ table is built.
		//objects->occ = &occ;
		//objects->EvaluateOcc(this->GetCurrentContainer().GetBlock().gt_ppa);

		for(uint32_t i = 0; i < objects->meta_container->size(); ++i){
			// Each entry evaluate occ if available.
			//entries[i].occ = objects->occ;
			//entries[i].EvaluateOcc();

			const uint32_t n_format_avail = entries[i].format_ids->size();
			if(n_format_avail > 0 && entries[i].is_loaded_gt){
				entries[i].gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				s.Update(entries[i], this->variant_container.GetAllocatedGenotypeMemory());
				//s.Update(entries[i]);
			}
		}

		delete [] entries;
		//objects->occ = nullptr;
		delete objects;
	}

	io::BasicBuffer json_buffer(250000);
	json_buffer += "{\"PSI\":{\n";
	for(uint32_t i = 0; i < this->GetGlobalHeader().GetNumberSamples(); ++i){
		s[i].LazyEvalute();
		if(i != 0) json_buffer += ",\n";
		s[i].ToJsonString(json_buffer, this->GetGlobalHeader().samples_[i]);
	}
	json_buffer += "\n}\n}\n";

	//std::cout << s.n_no_alts << "," << s.n_multi_allele << "," << s.n_multi_allele_snp << "," << s.n_biallelic << "," << s.n_singleton << std::endl;
	std::cout.write(json_buffer.data(), json_buffer.size());
	std::cout.flush();

	return true;
}

}
