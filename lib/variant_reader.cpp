#include "variant_reader.h"

#include "containers/variant_container.h"

#include "algorithm/compression/genotype_encoder.h"
#include "algorithm/permutation/genotype_sorter.h"

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
	char magic_string[TACHYON_MAGIC_HEADER_LENGTH];
	this->basic_reader.stream_.read(&magic_string[0], TACHYON_MAGIC_HEADER_LENGTH);
	if(strncmp(&magic_string[0], &TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH) != 0){
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

	io::BasicBuffer buf(100000);
	while(this->NextBlock()){
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		if(this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.GetBlock().gt_ppa);

		for(uint32_t i = 0; i < vc.size(); ++i){
			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if(this->GetBlockSettings().annotate_extra){
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}
			vc[i].ToVcfString(this->global_header, buf, this->GetBlockSettings().display_static, this->variant_container.GetAllocatedGenotypeMemory());
			std::cout.write(buf.data(), buf.size());
			buf.reset();
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

			if(this->GetBlockSettings().annotate_extra){
				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
			}

			vc[i].ToVcfString(this->global_header, buf, this->GetBlockSettings().display_static, this->variant_container.GetAllocatedGenotypeMemory());
			std::cout.write(buf.data(), buf.size());
			buf.reset();
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

	if(this->settings.annotate_genotypes)
		this->GetGlobalHeader().AddGenotypeAnnotationFields();

	// temp
	if(this->settings.group_file.size()){
		if(this->occ_table.ReadTable(this->settings.group_file, this->GetGlobalHeader(), '\t') == false){
			return(0);
		}

		this->global_header.AddGenotypeAnnotationFields(this->occ_table.row_names);
	}


	if(this->GetBlockSettings().show_vcf_header)
		this->global_header.PrintVcfHeader(std::cout);

	if(this->settings.use_htslib){
		if(this->interval_container.size()) return(this->OutputHtslibVcfSearch());
		else return(this->OutputHtslibVcfLinear());
	}

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

		if(this->GetBlockSettings().annotate_extra)
			this->occ_table.BuildTable(this->variant_container.GetBlock().gt_ppa);

		// Iterate over available records in this block.
		for(uint32_t i = 0; i < vc.size(); ++i){
			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if(this->GetBlockSettings().annotate_extra){
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			vc[i].OutputHtslibVcfInfo(rec, hdr, this->GetBlockSettings());
			vc[i].OutputHtslibVcfFormat(rec, hdr, this->GetBlockSettings(), this->GetCurrentContainer().GetAllocatedGenotypeMemory());
			vc[i].OutputHtslibVcfFilter(rec, hdr);


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

			if(this->GetBlockSettings().annotate_extra){
				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
			}

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			vc[i].OutputHtslibVcfInfo(rec, hdr, this->GetBlockSettings());
			vc[i].OutputHtslibVcfFormat(rec, hdr, this->GetBlockSettings(), this->GetCurrentContainer().GetAllocatedGenotypeMemory());
			vc[i].OutputHtslibVcfFilter(rec, hdr);

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

void VariantReader::UpdateHeaderView(void){
	io::VcfExtra e;
	e.key = "tachyon_viewVersion";
	e.value = TACHYON_PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  TACHYON_PROGRAM_NAME + '-' + TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	this->GetGlobalHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetGlobalHeader().extra_fields_.push_back(e);
	e.key = "tachyon_viewCommand";
	e.value = LITERAL_COMMAND_LINE;
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
			slave[i].global_header = &this->global_header;

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
	//const uint32_t n_threads = 1;
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
		slave_test[i].global_header = &this->global_header;

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
	s.Evaluate();
	s.ToJsonString(json_buffer, this->global_header.samples_);

	std::cout.write(json_buffer.data(), json_buffer.size());
	std::cout.flush();

	return true;
}

bool VariantReader::TempWrite(void){
	this->variant_container.AllocateGenotypeMemory();

	algorithm::GenotypeSorter  gts;
	gts.SetSamples(this->GetGlobalHeader().GetNumberSamples());
	algorithm::GenotypeEncoder gte(this->GetGlobalHeader().GetNumberSamples());

	//io::BasicBuffer buf(100000);
	while(this->NextBlock()){
		gts.reset();
		VariantContainer vc(this->GetCurrentContainer().GetBlock().header.n_variants);
		vc.Build(this->GetCurrentContainer().GetBlock(), this->global_header);

		containers::VariantBlock vblock;
		vblock.Allocate(100,100,100);
		vblock.resize(128000);

		if(this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.GetBlock().gt_ppa);

		for(uint32_t i = 0; i < vc.size(); ++i){
			if(this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if(this->GetBlockSettings().annotate_extra){
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}

			//vc[i].gt->Expand();
			vblock.AddMore(vc[i]);
			// Todo: Update controller (ref alt checks etc)
			vblock += vc[i];
		}

		vblock.header.controller.has_gt = false;
		vblock.header.controller.has_gt_permuted = false;

		if(gts.Build(vc.variants_, vc.n_variants_) == false){
			std::cerr << "failed to permute" << std::endl;
			return false;
		}

		if(gte.Encode(vc.variants_, vc.n_variants_, vc.block_, gts.permutation_array) == false){
			std::cerr << "failed to encode genotypes" << std::endl;
			return false;
		}


		//for(int i = 0; i < gts.permutation_array.n_s; ++i){
		//	std::cerr << "," << gts.permutation_array[i] << "=" << (*this->GetCurrentContainer().GetBlock().gt_ppa)[i];
		//}
		//std::cerr << std::endl;

		std::cerr << vblock.size() << ": " << vblock.header.minPosition << "-" << vblock.header.maxPosition << std::endl;
		vblock.UpdateContainers(this->global_header.GetNumberSamples());
		this->codec_manager.Compress(vblock,20,this->global_header.GetNumberSamples());

		std::cerr << "base: ";
		for(int i = 0; i < YON_BLK_N_STATIC; ++i){
			std::cerr << "," << YON_BLK_PRINT_NAMES[i] << ":" << vblock.base_containers[i].GetSizeUncompressed() << ":" << vblock.base_containers[i].header.HasMixedStride();
		}
		std::cerr << std::endl;

		std::cerr << "info: ";
		for(int i = 0; i < vblock.footer.n_info_streams; ++i){
			std::cerr << "," << i << ":" << vblock.info_containers[i].GetDataPrimitiveType() << ":" << vblock.info_containers[i].GetSizeUncompressed() << "->" << vblock.info_containers[i].GetSizeCompressed();
		}
		std::cerr << std::endl;

		std::cerr << "format: ";
		for(int i = 0; i < vblock.footer.n_format_streams; ++i){
			std::cerr << "," << i << ":" << vblock.format_containers[i].GetDataPrimitiveType() << ":" << vblock.format_containers[i].GetSizeUncompressed() << "->" << vblock.format_containers[i].GetSizeCompressed();
		}
		std::cerr << std::endl;

		vblock.Finalize();
		// After all compression and writing is finished the header
		// offsets are themselves compressed and stored in the block.
		vblock.PackFooter(); // Pack footer into buffer.
		this->codec_manager.zstd_codec.Compress(vblock.footer_support);
	}

	return 0;
}

}
