#include <fstream>
#include <regex>
#include <thread>

#include <openssl/evp.h>

#include "variant_writer.h"
#include "variant_importer.h"
#include "containers/checksum_container.h"
#include "algorithm/parallel/vcf_slaves.h"
#include "vcf_importer_slave.h"
#include "algorithm/timer.h"

namespace tachyon {

VariantImporterSettings::VariantImporterSettings() :
	verbose(true),
	permute_genotypes(true),
	encrypt_data(false),
	checkpoint_n_snps(1000),
	checkpoint_bases(10e6),
	n_threads(std::thread::hardware_concurrency()),
	compression_level(6),
	htslib_extra_threads(std::thread::hardware_concurrency() - 1 >= 0
	                     ? std::thread::hardware_concurrency() - 1
	                     : 0),
	info_end_key(-1),
	info_svlen_key(-1)
{
}

std::string VariantImporterSettings::GetInterpretedString(void) const {
	return(std::string("{\"input_file\":\"" + this->input_file +
	   "\",\"output_prefix\":\"" + this->output_prefix +
	   "\",\"checkpoint_snps\":" + std::to_string(this->checkpoint_n_snps) +
	   ",\"checkpoint_bases\":" + std::to_string(this->checkpoint_bases) +
	   ",\"compression_level\":" + std::to_string(this->compression_level) +
	   "}"
	));
}

/**<
 * PImpl implementation of private functionality of VariantImporter
 * class. Hides the private logic from the exposed public API and
 * speeds up compilation.
 */
class VariantImporter::VariantImporterImpl {
public:
	typedef VariantImporterImpl self_type;
	typedef io::VcfReader                 vcf_reader_type;
	typedef algorithm::CompressionManager compression_manager_type;
	typedef std::unordered_map<uint32_t, uint32_t> reorder_map_type;
	typedef std::unordered_map<uint64_t, uint32_t> hash_map_type;

	typedef VariantWriterInterface    writer_interface_type;
	typedef VariantWriterFile         writer_file_type;
	typedef VariantWriterStream       writer_stream_type;

public:
	~VariantImporterImpl() { delete this->writer; }
	bool Build(writer_interface_type* writer, settings_type& settings);

	void SetWriterTypeFile(void);
	void SetWriterTypeStream(void);

private:
	/**<
	 * If data has been encrypted then write out the provided encryption keychain
	 * using the destination writer.
	 * @param writer   Destination writer.
	 * @param keychain Source encryption keychain.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteKeychain(writer_interface_type* writer);

	/**<
	 * Write the appropriate yon archive header marker followed by the
	 * transmutation of the internal htslib Vcf header into a valid tachyon
	 * header and use the dst writer to write the output.
	 * @param writer Destination writer.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteYonHeader(writer_interface_type* writer);

	/**<
	 * Adds Extra fields to the global tachyon header describing when and how
	 * the target file was imported in the tachyon archive.
	 * @param header Global tachyon header.
	 */
	void UpdateHeaderImport(yon_vnt_hdr_t& header);

	/**<
	 * Converts a htslib-styled Vcf header into a Tachyon header.
	 * @param vcf_header Src reference of a Vcf header.
	 * @return           Returns a Tachyon yon_vnt_hdr_t.
	 */
	yon_vnt_hdr_t ConvertVcfHeader(const io::VcfHeader& vcf_header);

public:
	writer_interface_type* writer; // writer

	std::shared_ptr<settings_type> settings;
	compression_manager_type compression_manager;

	// Encryption keychain.
	Keychain keychain;

	// Map from BCF global FORMAT/INFO/FILTER IDX to local IDX such that
	// FORMAT maps to [0, f-1], and INFO maps to [0, i-1] and FILTER to
	// [0,l-1] and where f+i+l = n, where n is the total number of fields.
	//
	//                    Global    Local
	// std::unordered_map<uint32_t, uint32_t> filter_reorder_map_;
	reorder_map_type filter_reorder_map_;
	reorder_map_type info_reorder_map_;
	reorder_map_type format_reorder_map_;
	reorder_map_type contig_reorder_map_;

	std::unique_ptr<vcf_reader_type> vcf_reader_;
	yon_vnt_hdr_t yon_header_;

	hash_map_type block_hash_map;
};


VariantImporter::VariantImporter(void) :
	mImpl(new VariantImporter::VariantImporterImpl)
{

}

VariantImporter::VariantImporter(const settings_type& settings) :
	settings_(settings),
	mImpl(new VariantImporter::VariantImporterImpl)
{

}

VariantImporter::~VariantImporter() {
}

void VariantImporter::VariantImporterImpl::SetWriterTypeFile(void)  { this->writer = new writer_file_type;   }
void VariantImporter::VariantImporterImpl::SetWriterTypeStream(void) { this->writer = new writer_stream_type; }

bool VariantImporter::Build() {
	// Allocate a new writer.
	if (this->settings_.output_prefix.size() == 0 ||
	   (this->settings_.output_prefix.size() == 1 && this->settings_.output_prefix == "-"))
	{
		this->mImpl->writer = new VariantImporterImpl::writer_stream_type;
	}
	else this->mImpl->writer = new VariantImporterImpl::writer_file_type;

	// Open a file handle or standard out for writing.
	if (!this->mImpl->writer->open(this->settings_.output_prefix)) {
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	if (!this->mImpl->Build(this->mImpl->writer, this->settings_)) {
		std::cerr << utility::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
		return false;
	}
	return true;
}

bool VariantImporter::VariantImporterImpl::Build(writer_interface_type* writer, settings_type& settings) {
	if (writer == nullptr)
		return false;

	this->settings = std::shared_ptr<settings_type>(&settings);

	// Retrieve a unique VcfReader.
	this->vcf_reader_ = io::VcfReader::FromFile(this->settings->input_file, this->settings->htslib_extra_threads);
	if (this->vcf_reader_ == nullptr)
		return false;


	for (uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i) {
		if (this->vcf_reader_->vcf_header_.contigs_[i].n_bases == 0) {
			std::cerr << utility::timestamp("NOTICE") << "No length declared for contig (. " << this->vcf_reader_->vcf_header_.contigs_[i].name << "). Setting to INT32_MAX." << std::endl;
			this->vcf_reader_->vcf_header_.contigs_[i].n_bases = std::numeric_limits<int32_t>::max();
		}
	}

	// Remap the global IDX fields in Vcf to the appropriate incremental order.
	// This is useful in the situations when fields have been removed or added
	// to the Vcf header section without reformatting the file.
	for (uint32_t i = 0; i < this->vcf_reader_->vcf_header_.contigs_.size(); ++i)
		this->contig_reorder_map_[this->vcf_reader_->vcf_header_.contigs_[i].idx] = i;

	for (uint32_t i = 0; i < this->vcf_reader_->vcf_header_.info_fields_.size(); ++i)
		this->info_reorder_map_[this->vcf_reader_->vcf_header_.info_fields_[i].idx] = i;

	for (uint32_t i = 0; i < this->vcf_reader_->vcf_header_.format_fields_.size(); ++i)
		this->format_reorder_map_[this->vcf_reader_->vcf_header_.format_fields_[i].idx] = i;

	for (uint32_t i = 0; i < this->vcf_reader_->vcf_header_.filter_fields_.size(); ++i)
		this->filter_reorder_map_[this->vcf_reader_->vcf_header_.filter_fields_[i].idx] = i;

	// Predicate of a search for "GT" FORMAT field in the Vcf header.
	bool GT_available = (this->vcf_reader_->vcf_header_.GetFormat("GT") != nullptr);

	// Predicate of a search for "END" INFO field in the Vcf header.
	VcfInfo* vcf_info_end = this->vcf_reader_->vcf_header_.GetInfo("END");
	if (vcf_info_end != nullptr)
		this->settings->info_end_key = vcf_info_end->idx;

	// Setup the checksums container.
	algorithm::VariantDigestManager checksums(YON_BLK_N_STATIC   + 1, // Add one for global checksum.
			this->vcf_reader_->vcf_header_.info_fields_.size()   + 1,
			this->vcf_reader_->vcf_header_.format_fields_.size() + 1);

	// The index needs to know how many contigs that's described in the
	// Vcf header and their lengths in base-pairs. This information is
	// needed to construct the linear and quad-tree index most appropriate
	// for the data lengths.
	writer->index.Setup(this->vcf_reader_->vcf_header_.contigs_);

	// Write out a fresh Tachyon header with the data from the Vcf header. As
	// this data will not be modified during the import stage it is safe to
	// write out now.
	this->WriteYonHeader(writer);

	yon_producer_vcfc  producer(*this->settings, this->vcf_reader_, this->settings->n_threads);
	yon_consumer_vcfc* consumers = new yon_consumer_vcfc[this->settings->n_threads];
	yon_writer_sync    write;
	write.writer = writer;

	write.stats_basic.Allocate(YON_BLK_N_STATIC);
	write.stats_info.Allocate(this->vcf_reader_->vcf_header_.info_fields_.size());
	write.stats_format.Allocate(this->vcf_reader_->vcf_header_.format_fields_.size());

	// Start progress timer.
	algorithm::Timer timer;
	timer.Start();

	// Start producer thread with multiple htslib threads.
	producer.Start();
	// Setup shared pointers for slaves.
	std::shared_ptr<io::VcfHeader> s_vcf_hdr(&this->vcf_reader_->vcf_header_);
	std::shared_ptr<VariantImporterSettings> s_settings(this->settings);
	std::shared_ptr<std::atomic<bool>> s_data_avail(&producer.data_available);
	std::shared_ptr<yon_pool_vcfc> s_data_pool(&producer.data_pool);
	std::shared_ptr<yon_writer_sync> s_writer(&write);
	std::shared_ptr<Keychain> s_keychain(&this->keychain);

	// Setup and start slaves.
	for (uint32_t i = 0; i < this->settings->n_threads; ++i) {
		consumers[i].thread_id = i;
		consumers[i].data_available = s_data_avail;
		consumers[i].data_pool = s_data_pool;
		consumers[i].global_header = s_vcf_hdr;
		consumers[i].poolw = s_writer;
		consumers[i].importer.GT_available_       = GT_available;
		consumers[i].importer.SetVcfHeader(s_vcf_hdr);
		consumers[i].importer.settings_           = s_settings;
		consumers[i].importer.keychain            = s_keychain;
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

	// Join consumer and producer threads.
	for (uint32_t i = 0; i < this->settings->n_threads; ++i) consumers[i].thread_.join();
	producer.all_finished = true;
	producer.thread_.join();

	// Reduce functions.
	for (uint32_t i = 1; i < this->settings->n_threads; ++i) consumers[0] += consumers[i];
	write.stats_basic  += consumers[0].importer.stats_basic;
	write.stats_info   += consumers[0].importer.stats_info;
	write.stats_format += consumers[0].importer.stats_format;

	writer->index += consumers[0].importer.index;
	//this->writer->index.Print(std::cerr);

	// Finalize writing procedure.
	writer->close();
	this->WriteKeychain(writer);

	uint64_t b_uncompressed = 0;
	if (settings.output_prefix != "-") {
		std::cout << "Field\tType\tCompressed\tUncompressed\tStrideCompressed\tStrideUncompressed\tFold\tBinaryVcf\tFold-Bcf-Yon" << std::endl;
		for (int i = 0; i < write.stats_basic.size(); ++i) {
			std::cout << YON_BLK_PRINT_NAMES[i] << "\tNA\t" << write.stats_basic.at(i) << std::endl;
			b_uncompressed += write.stats_basic[i].cost_uncompressed + write.stats_basic[i].cost_strides;
		}

		for (int i = 0; i < write.stats_info.size(); ++i) {
			std::cout << "INFO-" << this->yon_header_.info_fields_[i].id << "\t" << this->yon_header_.info_fields_[i].type << "\t" << write.stats_info.at(i) << std::endl;
			b_uncompressed += write.stats_info[i].cost_uncompressed + write.stats_info[i].cost_strides;
		}

		for (int i = 0; i < write.stats_format.size(); ++i) {
			std::cout << "FORMAT-" << this->yon_header_.format_fields_[i].id << "\t" << this->yon_header_.format_fields_[i].type << "\t" << write.stats_format.at(i) << std::endl;
			b_uncompressed += write.stats_format[i].cost_uncompressed + write.stats_format[i].cost_strides;
		}
	}

	std::cerr << utility::timestamp("PROGRESS") << "Processed " << utility::ToPrettyDiskString(consumers[0].b_indiv + consumers[0].b_shared) << " of htslib bcf1_t records -> " << utility::ToPrettyDiskString(b_uncompressed) << " (" << (double)(consumers[0].b_indiv + consumers[0].b_shared)/b_uncompressed << ")" << std::endl;
	std::cerr << utility::timestamp("PROGRESS") << "Wrote: " << utility::ToPrettyString(consumers[0].poolw->n_written_rcds) << " variants to " << utility::ToPrettyString(consumers->poolw->writer->n_blocks_written) << " blocks in " << utility::ToPrettyDiskString((uint64_t)writer->stream->tellp()) << "(" << (double)(consumers[0].b_indiv + consumers[0].b_shared)/((uint64_t)writer->stream->tellp()) << ")" << std::endl;
	std::cerr << utility::timestamp("PROGRESS") << "All done (" << timer.ElapsedString() << ")" << std::endl;

	//delete [] consumers;

	// All done
	return(true);
}

bool VariantImporter::VariantImporterImpl::WriteKeychain(writer_interface_type* writer) {
	// Write encryption keychain.
	if (this->settings->encrypt_data) {
		if (this->settings->output_prefix.size()) {
			std::ofstream writer_keychain;
			writer_file_type* wstats = reinterpret_cast<writer_file_type*>(writer);
			writer_keychain.open(wstats->basePath + wstats->baseName + ".kyon", std::ios::out);
			if (!SILENT)
				std::cerr << utility::timestamp("LOG") << "Writing encryption keychain to: " << (wstats->basePath + wstats->baseName) << ".kyon" << std::endl;

			if (writer_keychain.good()) {
				//writer_keychain.write(keybuffer.data(), keybuffer.size());
				writer_keychain << keychain;
				writer_keychain.flush();
			}
			const uint32_t keychain_size = writer_keychain.tellp();
			writer_keychain.close();
			if (!SILENT)
				std::cerr << utility::timestamp("LOG") << "Wrote keychain with " << utility::ToPrettyString(keychain.size()) << " keys to " << utility::ToPrettyDiskString(keychain_size) << "..." << std::endl;
		}
	}
	return true;
}

bool VariantImporter::VariantImporterImpl::WriteYonHeader(writer_interface_type* writer) {
	if (writer == nullptr)
		return false;

	// Write basic header prefix.
	writer->stream->write(&TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH); // Todo: fix

	// Transmute a htslib-styled vcf header into a tachyon
	// header.
	this->yon_header_ = this->ConvertVcfHeader(this->vcf_reader_->vcf_header_);
	//this->yon_header_ = yon_vnt_hdr_t(this->vcf_reader_->vcf_header_);
	// Update the extra provenance fields in the new header.
	this->UpdateHeaderImport(this->yon_header_);

	// Pack header into a byte-stream, compress it, and write
	// it out.
	yon_buffer_t temp(500000);
	yon_buffer_t temp_cmp(temp);
	temp << this->yon_header_;
	this->compression_manager.zstd_codec.Compress(temp, temp_cmp, 20);
	uint32_t l_data   = temp.size();
	uint32_t l_c_data = temp_cmp.size();
	utility::SerializePrimitive(l_data,   *writer->stream);
	utility::SerializePrimitive(l_c_data, *writer->stream);
	writer->stream->write(temp_cmp.data(), l_c_data);
	return(writer->stream->good());
}

void VariantImporter::VariantImporterImpl::UpdateHeaderImport(yon_vnt_hdr_t& header) {
	VcfExtra e;
	e.key = "tachyon_importVersion";
	e.value = tachyon::TACHYON_PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  tachyon::TACHYON_PROGRAM_NAME + '-' + tachyon::TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importCommand";
	e.value = tachyon::LITERAL_COMMAND_LINE;
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importSettings";
	e.value = this->settings->GetInterpretedString();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
}

yon_vnt_hdr_t VariantImporter::VariantImporterImpl::ConvertVcfHeader(const io::VcfHeader& vcf_header) {
	yon_vnt_hdr_t header;
	header.fileformat_string_ = vcf_header.fileformat_string_;
	header.literals_ = vcf_header.literals_;
	header.samples_ = vcf_header.samples_;
	header.filter_fields_ = vcf_header.filter_fields_;
	header.structured_extra_fields_ = vcf_header.structured_extra_fields_;
	header.extra_fields_ = vcf_header.extra_fields_;

	header.BuildMaps();
	header.BuildReverseMaps();

	header.contigs_.resize(vcf_header.contigs_.size());
	for (uint32_t i = 0; i < vcf_header.contigs_.size(); ++i)
		header.contigs_[i] = vcf_header.contigs_[i];

	header.info_fields_.resize(vcf_header.info_fields_.size());
	for (uint32_t i = 0; i < vcf_header.info_fields_.size(); ++i)
		header.info_fields_[i] = vcf_header.info_fields_[i];

	header.format_fields_.resize(vcf_header.format_fields_.size());
	for (uint32_t i = 0; i < vcf_header.format_fields_.size(); ++i)
		header.format_fields_[i] = vcf_header.format_fields_[i];

	header.RecodeIndices();

	return(header);
}

} /* namespace Tachyon */
