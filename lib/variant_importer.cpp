#include <fstream>
#include <regex>

#include "variant_importer.h"
#include "containers/checksum_container.h"

#include "algorithm/parallel/vcf_slaves.h"
#include "vcf_importer_slave.h"

namespace tachyon {

VariantImporter::VariantImporter(const settings_type& settings) :
		settings_(settings),
		writer(nullptr)
{

}

VariantImporter::~VariantImporter(){
	delete this->writer;
}

void VariantImporter::clear(){
	//
}

bool VariantImporter::Build(){
	//if(!this->BuildVCF()){
	if(!this->BuildParallel()){
		std::cerr << utility::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
		return false;
	}
	return true;
}

/*
bool VariantImporter::BuildVCF(void){
	// Retrieve a unique VcfReader.
	this->vcf_reader_ = io::VcfReader::FromFile(this->settings_.input_file);
	if(this->vcf_reader_ == nullptr){
		return false;
	}

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i){
		if(this->vcf_reader_->vcf_header_.contigs_[i].n_bases == 0){
			std::cerr << utility::timestamp("NOTICE") << "No length declared for contig. Setting to INT32_MAX." << std::endl;
			this->vcf_reader_->vcf_header_.contigs_[i].n_bases = std::numeric_limits<int32_t>::max();
		}
	}

	// Remap the global IDX fields in Vcf to the appropriate incremental order.
	// This is useful in the situations when fields have been removed or added
	// to the Vcf header section without reformatting the file.
	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i)
		this->contig_reorder_map_[this->vcf_reader_->vcf_header_.contigs_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.info_fields_.size(); ++i)
		this->info_reorder_map_[this->vcf_reader_->vcf_header_.info_fields_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.format_fields_.size(); ++i)
		this->format_reorder_map_[this->vcf_reader_->vcf_header_.format_fields_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.filter_fields_.size(); ++i)
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
	const uint32_t resize_to = this->settings_.checkpoint_n_snps * sizeof(uint32_t) * 2; // small initial allocation
	this->block.resize(resize_to);

	// Start porgress timer.
	algorithm::Timer timer; timer.Start();

	// Iterate over all available variants in the file or until encountering
	// an error.
	while(true){
		// Retrieve bcf1_t records using htslib and lazy evaluate them. Stop
		// after retrieving a set number of variants or if the interval between
		// the smallest and largest variant exceeds some distance in base pairs.
		if(this->vcf_container_.GetVariants(this->settings_.checkpoint_n_snps,
		                                    this->settings_.checkpoint_bases,
		                                    this->vcf_reader_) == false)
		{
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
*/

bool VariantImporter::BuildParallel(void){
	// Retrieve a unique VcfReader.
	this->vcf_reader_ = io::VcfReader::FromFile(this->settings_.input_file, this->settings_.htslib_extra_threads);
	if(this->vcf_reader_ == nullptr){
		return false;
	}

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i){
		if(this->vcf_reader_->vcf_header_.contigs_[i].n_bases == 0){
			std::cerr << utility::timestamp("NOTICE") << "No length declared for contig (. " << this->vcf_reader_->vcf_header_.contigs_[i].name << "). Setting to INT32_MAX." << std::endl;
			this->vcf_reader_->vcf_header_.contigs_[i].n_bases = std::numeric_limits<int32_t>::max();
		}
	}

	// Remap the global IDX fields in Vcf to the appropriate incremental order.
	// This is useful in the situations when fields have been removed or added
	// to the Vcf header section without reformatting the file.
	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i)
		this->contig_reorder_map_[this->vcf_reader_->vcf_header_.contigs_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.info_fields_.size(); ++i)
		this->info_reorder_map_[this->vcf_reader_->vcf_header_.info_fields_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.format_fields_.size(); ++i)
		this->format_reorder_map_[this->vcf_reader_->vcf_header_.format_fields_[i].idx] = i;

	for(uint32_t i = 0; i < this->vcf_reader_->vcf_header_.filter_fields_.size(); ++i)
		this->filter_reorder_map_[this->vcf_reader_->vcf_header_.filter_fields_[i].idx] = i;

	// Predicate of a search for "GT" FORMAT field in the Vcf header.
	bool GT_available = (this->vcf_reader_->vcf_header_.GetFormat("GT") != nullptr);

	// Predicate of a search for "END" INFO field in the Vcf header.
	io::VcfInfo* vcf_info_end = this->vcf_reader_->vcf_header_.GetInfo("END");
	if(vcf_info_end != nullptr)
		this->settings_.info_end_key = vcf_info_end->idx;

	// Setup the checksums container.
	algorithm::VariantDigestManager checksums(YON_BLK_N_STATIC   + 1, // Add one for global checksum.
			this->vcf_reader_->vcf_header_.info_fields_.size()   + 1,
			this->vcf_reader_->vcf_header_.format_fields_.size() + 1);

	// Setup the encryption container.
	encryption::EncryptionDecorator encryption_manager;
	encryption::Keychain<> keychain;

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

	// The index needs to know how many contigs that's described in the
	// Vcf header and their lengths in base-pairs. This information is
	// needed to construct the linear and quad-tree index most appropriate
	// for the data lengths.
	this->writer->index.Setup(this->vcf_reader_->vcf_header_.contigs_);

	// Write out a fresh Tachyon header with the data from the Vcf header. As
	// this data will not be modified during the import stage it is safe to
	// write out now.
	this->WriteYonHeader();

	yon_producer_vcfc  producer(this->settings_, this->vcf_reader_, this->settings_.n_threads);
	yon_consumer_vcfc* consumers = new yon_consumer_vcfc[this->settings_.n_threads];
	yon_writer_sync    write;
	write.writer = this->writer;

	write.stats_basic.Allocate(YON_BLK_N_STATIC);
	write.stats_info.Allocate(this->vcf_reader_->vcf_header_.info_fields_.size());
	write.stats_format.Allocate(this->vcf_reader_->vcf_header_.format_fields_.size());

	// Start progress timer.
	algorithm::Timer timer; timer.Start();
	producer.Start();
	for(uint32_t i = 0; i < this->settings_.n_threads; ++i){
		consumers[i].thread_id = i;
		consumers[i].data_available = &producer.data_available;
		consumers[i].data_pool = &producer.data_pool;
		consumers[i].global_header = &this->vcf_reader_->vcf_header_;
		consumers[i].poolw = &write;
		consumers[i].importer.GT_available_       = GT_available;
		consumers[i].importer.SetVcfHeader(&this->vcf_reader_->vcf_header_);
		consumers[i].importer.settings_           = this->settings_;
		consumers[i].importer.info_reorder_map_   = this->info_reorder_map_;
		consumers[i].importer.format_reorder_map_ = this->format_reorder_map_;
		consumers[i].importer.filter_reorder_map_ = this->filter_reorder_map_;
		consumers[i].importer.contig_reorder_map_ = this->contig_reorder_map_;

		consumers[i].importer.stats_basic.Allocate(YON_BLK_N_STATIC);
		consumers[i].importer.stats_info.Allocate(this->vcf_reader_->vcf_header_.info_fields_.size());
		consumers[i].importer.stats_format.Allocate(this->vcf_reader_->vcf_header_.format_fields_.size());

		// The index needs to know how many contigs that's described in the
		// Vcf header and their lenghts. This information is needed to construct
		// the linear and quad-tree index most appropriate for the data.
		consumers[i].importer.index.Setup(this->vcf_reader_->vcf_header_.contigs_);
		consumers[i].Start();
	}

	for(uint32_t i = 0; i < this->settings_.n_threads; ++i) consumers[i].thread_.join();
	producer.all_finished = true;
	producer.thread_.join();
	for(uint32_t i = 1; i < this->settings_.n_threads; ++i) consumers[0] += consumers[i];
	write.stats_basic += consumers[0].importer.stats_basic;
	write.stats_info += consumers[0].importer.stats_info;
	write.stats_format += consumers[0].importer.stats_format;

	//std::cerr << "blocks processed: " << consumers[0].importer.n_blocks_processed << std::endl;
	//std::cerr << producer.n_rcds_loaded << "==" << consumers[0].n_rcds_processed << std::endl;
	//std::cerr << "processed: " << consumers[0].b_indiv << "b and " << consumers[0].b_shared << std::endl;
	std::cerr << utility::timestamp("PROGRESS") << "Processed " << utility::toPrettyDiskString(consumers[0].b_indiv + consumers[0].b_shared) << " of htslib bcf1_t records." << std::endl;
	//std::cerr << "wrote: " << consumers[0].poolw->n_written_rcds << "rcds to " << consumers->poolw->writer->n_blocks_written << " writer says " << consumers->poolw->writer->n_variants_written << std::endl;


	std::cerr << "Field\tCompressed\tUncompressed\tStrideCompressed\tStrideUncompressed\tFold\tBinaryVcf\tFold-Bcf-Yon" << std::endl;
	for(int i = 0; i < write.stats_basic.size(); ++i)
		std::cerr << YON_BLK_PRINT_NAMES[i] << "\t" << write.stats_basic.at(i) << std::endl;

	for(int i = 0; i < write.stats_info.size(); ++i)
		std::cerr << "INFO-" << this->yon_header_.info_fields_[i].id << "\t" << write.stats_info.at(i) << std::endl;

	for(int i = 0; i < write.stats_format.size(); ++i)
		std::cerr << "FORMAT-" << this->yon_header_.format_fields_[i].id << "\t" << write.stats_format.at(i) << std::endl;


	this->writer->index += consumers[0].importer.index;
	//this->writer->index.Print(std::cerr);

	// Finalize writing procedure.
	write.WriteFinal(checksums);
	//write.WriteKeychain(keychain);

	std::cerr << utility::timestamp("PROGRESS") << "All done (" << timer.ElapsedString() << ")" << std::endl;

	// All done
	return(true);
}

bool VariantImporter::WriteKeychain(const encryption::Keychain<>& keychain){
	// Write encryption keychain.
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
			const uint32_t keychain_size = writer_keychain.tellp();
			writer_keychain.close();
			if(!SILENT)
				std::cerr << utility::timestamp("LOG") << "Wrote keychain with " << utility::ToPrettyString(keychain.size()) << " keys to " << utility::toPrettyDiskString(keychain_size) << "..." << std::endl;
		}
	}
	return true;
}

bool VariantImporter::WriteYonHeader(){
	// Write basic header prefix.
	this->writer->stream->write(&constants::FILE_HEADER[0], constants::FILE_HEADER_LENGTH); // Todo: fix

	// Transmute a htslib-styled vcf header into a tachyon
	// header.
	this->yon_header_ = VariantHeader(this->vcf_reader_->vcf_header_);
	// Update the extra provenance fields in the new header.
	this->UpdateHeaderImport(this->yon_header_);

	// Pack header into a byte-stream, compress it, and write
	// it out.
	io::BasicBuffer temp(500000);
	io::BasicBuffer temp_cmp(temp);
	temp << this->yon_header_;
	this->compression_manager.zstd_codec.Compress(temp, temp_cmp, 20);
	uint32_t l_data   = temp.size();
	uint32_t l_c_data = temp_cmp.size();
	utility::SerializePrimitive(l_data,   *this->writer->stream);
	utility::SerializePrimitive(l_c_data, *this->writer->stream);
	this->writer->stream->write(temp_cmp.data(), l_c_data);
	return(this->writer->stream->good());
}

void VariantImporter::UpdateHeaderImport(VariantHeader& header){
	io::VcfExtra e;
	e.key = "tachyon_importVersion";
	e.value = tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importCommand";
	e.value = tachyon::constants::LITERAL_COMMAND_LINE;
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importSettings";
	e.value = this->settings_.GetInterpretedString();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
}

} /* namespace Tachyon */
