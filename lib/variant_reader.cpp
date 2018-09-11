#include "variant_reader.h"

#include "containers/variant_container.h"

namespace tachyon {

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

	// Validate the global footer.
	if(this->global_footer.Validate() == false){
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

bool VariantReader::NextBlockRaw(void){
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

uint64_t VariantReader::OutputVcfLinear(void){
	this->variant_container.AllocateGenotypeMemory();
	// temp
	//if(this->occ_table.ReadTable("/media/mdrk/NVMe/1kgp3/populations/integrated_call_samples_v3.20130502.ALL.panel", this->GetGlobalHeader(), '\t') == false){
	//	return(0);
	//}

	io::BasicBuffer buf(100000);
	while(this->NextBlock()){
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		for(uint32_t i = 0; i < vc.size(); ++i){
			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			for(int i = 0; i < vc.n_variants_; ++i){
				vc[i].Print(this->global_header, buf, this->GetBlockSettings().display_static, this->variant_container.GetAllocatedGenotypeMemory());
				std::cout.write(buf.data(), buf.size());
				buf.reset();
			}
		}
	}

	return 0;
}

uint64_t VariantReader::OutputVcfSearch(void){
	this->variant_container.AllocateGenotypeMemory();

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;
	io::BasicBuffer buf(100000);

	for(uint32_t i = 0; i < this->interval_container.GetBlockList().size(); ++i){
		this->GetBlock(this->interval_container.GetBlockList()[i]);

		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);


		for(uint32_t i = 0; i < vc.size(); ++i){
			if((this->*filter_intervals)(vc[i]) == false)
				continue;

			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			for(int i = 0; i < vc.n_variants_; ++i){
				vc[i].Print(this->global_header, buf, this->GetBlockSettings().display_static, this->variant_container.GetAllocatedGenotypeMemory());
				std::cout.write(buf.data(), buf.size());
				buf.reset();
			}
		}
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
	hts_stream_type[0] = 'w';
	hts_stream_type[1] = this->settings.output_type;
	htsFile* fp = hts_open(this->settings.output.c_str(), hts_stream_type);

	// Add extra htslib compression threads if desired.
	int n_extra_threads = std::thread::hardware_concurrency();
	if(n_extra_threads){
		int ret = hts_set_threads(fp, n_extra_threads);
		if(ret < 0){
			std::cerr << "failed to open multiple handles" << std::endl;
			return 0;
		}
	}

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
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		// Iterate over available records in this block.
		for(uint32_t i = 0; i < vc.size(); ++i){
			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			this->OutputHtslibVcfInfo(rec, hdr, vc[i]);
			this->OutputHtslibVcfFormat(rec, hdr, vc[i]);
			this->OutputHtslibVcfFilter(rec, hdr, vc[i]);

			if ( bcf_write1(fp, hdr, rec) != 0 ){
				std::cerr << "Failed to write record to " << this->settings.output;
				exit(1);
			}

			bcf_clear1(rec);
		}
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

	// Add extra htslib compression threads if desired.
	int n_extra_threads = std::thread::hardware_concurrency();
	if(n_extra_threads){
		int ret = hts_set_threads(fp, n_extra_threads);
		if(ret < 0){
			std::cerr << "failed to open multiple handles" << std::endl;
			return 0;
		}
	}

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
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		// Iterate over available records in this block.
		for(uint32_t i = 0; i < vc.size(); ++i){
			if((this->*filter_intervals)(vc[i]) == false)
				continue;

			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			this->OutputHtslibVcfInfo(rec, hdr, vc[i]);
			this->OutputHtslibVcfFormat(rec, hdr, vc[i]);
			this->OutputHtslibVcfFilter(rec, hdr, vc[i]);

			if ( bcf_write1(fp, hdr, rec) != 0 ){
				std::cerr << "Failed to write record to " << this->settings.output;
				exit(1);
			}

			bcf_clear1(rec);
		}
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

void VariantReader::OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr, yon1_vnt_t& entry) const{
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

void VariantReader::OutputHtslibVcfFormat(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_vnt_t& entry) const{
	if(entry.n_fmt){
		const uint32_t n_format_avail = entry.fmt_ids->size();
		if(n_format_avail){
			// Case when the only available FORMAT field is the GT field.
			if(n_format_avail == 1 && entry.is_loaded_gt &&
			   entry.controller.gt_available &&
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
			        entry.controller.gt_available &&
			        (this->GetBlockSettings().display_static & YON_BLK_BV_GT))
			{
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				entry.gt->UpdateHtslibGenotypes(rec, hdr);
				for(uint32_t g = 1; g < n_format_avail; ++g){
					if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, entry.fmt_hdr[g]->id);
					else if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, entry.fmt_hdr[g]->id);
					else if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, entry.fmt_hdr[g]->id);
				}
				entry.gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				for(uint32_t g = 0; g < n_format_avail; ++g){
					if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, entry.fmt_hdr[g]->id);
					else if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, entry.fmt_hdr[g]->id);
					else if(entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || entry.fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
						entry.fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, entry.fmt_hdr[g]->id);
				}
			}
		}
	}
}

void VariantReader::OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_vnt_t& entry) const{
	if(entry.n_flt){
		for(uint32_t k = 0; k < entry.flt_ids->size(); ++k){
			int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, entry.flt_hdr[k]->id.data());
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

bool VariantReader::Benchmark(const uint32_t threads){
	std::cerr << utility::timestamp("LOG") << "Starting benchmark with " << threads << " threads..." << std::endl;
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::LoadData);
	this->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::UncompressData);
	this->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::EvaluateData);
	this->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::EvaluateRecords);
	this->basic_reader.stream_.close();
	return(true);
}

bool VariantReader::BenchmarkWrapper(const uint32_t threads, bool(VariantSlavePerformance::*func)(containers::VariantBlock*&)){
	if(!this->open(settings.input)){
			std::cerr << "failed to open" << std::endl;
			return 1;
		}

		this->GetBlockSettings().LoadAll(true);

		// Begin
		algorithm::Timer timer; timer.Start();
		const uint32_t n_threads = threads;
		yon_producer_vblock<VariantReader> prd(n_threads);
		prd.Setup(&VariantReader::NextBlockRaw, *this, this->variant_container.GetBlock());
		prd.Start();

		yon_consumer_vblock<VariantSlavePerformance>* csm = new yon_consumer_vblock<VariantSlavePerformance>[n_threads];
		VariantSlavePerformance* slave = new VariantSlavePerformance[n_threads];

		for(int i = 0; i < n_threads; ++i){
			csm[i].thread_id      = i;
			csm[i].data_available = &prd.data_available;
			csm[i].data_pool      = &prd.data_pool;

			slave[i].settings = this->GetBlockSettings();
			slave[i].vc << this->global_header;
			slave[i].vc.AllocateGenotypeMemory();

			csm[i].Start(func, slave[i]);
		}
		// Join consumer and producer threads.
		for(uint32_t i = 0; i < n_threads; ++i) csm[i].thread_.join();

		prd.all_finished = true;
		prd.thread_.join();

		for(uint32_t i = 1; i < n_threads; ++i) slave[0] += slave[i];

		delete [] slave;
		delete [] csm;

		const double time_elapsed = timer.Elapsed().count();
		std::cerr << utility::timestamp("LOG") << time_elapsed << "\t" << slave[0].data_loaded << "\t" << (double)slave[0].data_loaded/time_elapsed/1e6 << "\t" << slave[0].data_uncompressed << "\t" << (double)slave[0].data_uncompressed/time_elapsed/1e6 << std::endl;
		return true;
}

bool VariantReader::Stats(void){
	const uint32_t n_threads = std::thread::hardware_concurrency();
	yon_producer_vblock<VariantReader> prd(n_threads);
	prd.Setup(&VariantReader::NextBlockRaw, *this, this->variant_container.GetBlock());
	prd.Start();

	yon_consumer_vblock<VariantSlaveTsTv>* csm = new yon_consumer_vblock<VariantSlaveTsTv>[n_threads];

	VariantSlaveTsTv* slave_test = new VariantSlaveTsTv[n_threads];

 	for(int i = 0; i < n_threads; ++i){
		csm[i].thread_id      = i;
		csm[i].data_available = &prd.data_available;
		csm[i].data_pool      = &prd.data_pool;

		slave_test[i].settings = this->GetBlockSettings();
		slave_test[i].s.SetSize(this->global_header.GetNumberSamples());
		slave_test[i].s_local.SetSize(this->global_header.GetNumberSamples());
		slave_test[i].vc << this->global_header;
		slave_test[i].vc.AllocateGenotypeMemory();

		csm[i].Start(&VariantSlaveTsTv::GatherGenotypeStatistics, slave_test[i]);
	}
	// Join consumer and producer threads.
	for(uint32_t i = 0; i < n_threads; ++i) csm[i].thread_.join();

	prd.all_finished = true;
	prd.thread_.join();

	yon_stats_tstv s(this->GetGlobalHeader().GetNumberSamples());
	for(uint32_t i = 0; i < n_threads; ++i) s += slave_test[i].s;

	delete [] slave_test;
	delete [] csm;

	io::BasicBuffer json_buffer(250000);
	s.ToJsonString(json_buffer, this->global_header.samples_);

	std::cout.write(json_buffer.data(), json_buffer.size());
	std::cout.flush();

	return true;
}

bool VariantReader::TempWrite(void){
	io::BasicBuffer buf(256000);
	this->variant_container.AllocateGenotypeMemory();

	while(this->NextBlock()){
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		for(int i = 0; i < vc.n_variants_; ++i){
			vc[i].Print(this->global_header, buf, this->GetBlockSettings().display_static, this->variant_container.GetAllocatedGenotypeMemory());
			std::cout.write(buf.data(), buf.size());
			buf.reset();
		}

		//delete objects;
	}

	return 0;
}

}
