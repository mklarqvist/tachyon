#include <variant_container.h>
#include <openssl/evp.h>

#include "variant_reader.h"
#include "variant_writer.h"

#include "algorithm/compression/genotype_encoder.h"
#include "algorithm/permutation/genotype_sorter.h"
#include "algorithm/timer.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "containers/interval_container.h"
#include "io/basic_reader.h"

#include "algorithm/parallel/vcf_slaves.h"
#include "algorithm/parallel/variant_slaves.h"
#include "algorithm/parallel/variant_base_slave.h"

namespace tachyon {

VariantReaderSettings::VariantReaderSettings() :
	drop_format(false), header_only(false), show_header(true),
	annotate_genotypes(false), use_htslib(false),
	output("-"), output_type('v'),
	permute_genotypes(true), encrypt_data(true),
	checkpoint_n_snps(500), checkpoint_bases(5000000),
	n_threads(std::thread::hardware_concurrency()), compression_level(6)
{}

std::string VariantReaderSettings::GetSettingsString(void) const {
	return(std::string(
	"{"
	"\"input\":" + (this->input.length() ? "\"" + this->input + "\"" : "null") +
	",\"output\":" + (this->output.length() ? "\"" + this->output + "\"" : "null") +
	",\"keychain_file\":" + (this->keychain_file.length() ? "\"" + this->keychain_file + "\"" : "null") +
	",\"annotate_genotypes\":" + (this->annotate_genotypes ? "true" : "false") +
	",\"drop_format\":" + (this->drop_format ? "true" : "false") +
	"}; timestamp=" + tachyon::utility::datetime()
	));
}

class VariantReader::VariantReaderImpl {
public:
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef algorithm::VariantDigestManager        checksum_type;
	typedef containers::IntervalContainer          interval_container_type;
	typedef algorithm::Interval<uint32_t, int64_t> interval_type;
	typedef io::BasicReader                        basic_reader_type;

public:
	VariantReaderImpl() = default;
	VariantReaderImpl(const std::string& filename) : basic_reader(filename) {}

	/**<
	* Converts this header object into a hts_vcf_header object from the
	* internally stored literal string. This object is required for
	* writing out VCF/BCF files.
	* @return
	*/
	bcf_hdr_t* ConvertVcfHeaderLiterals(const yon_vnt_hdr_t& hdr, const bool add_format);
	bcf_hdr_t* ConvertVcfHeader(const yon_vnt_hdr_t& hdr, const bool add_format);

public:
	basic_reader_type       basic_reader;
	checksum_type           checksums;
	codec_manager_type      codec_manager;
	interval_container_type interval_container;

};

VariantReader::VariantReader() :
	mImpl(new VariantReader::VariantReaderImpl),
	b_data_start(0),
	gt_exp(nullptr)
{}

VariantReader::VariantReader(const std::string& filename) :
	mImpl(new VariantReader::VariantReaderImpl(filename)),
	b_data_start(0),
	gt_exp(nullptr)
{

}

VariantReader::~VariantReader() {
	delete [] gt_exp;
}

bool VariantReader::open(void) {
	if (this->settings.input.size() == 0) {
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	if (this->mImpl->basic_reader.open() == false) {
		std::cerr << "Failed to open" << std::endl;
		return false;
	}

	if (this->mImpl->basic_reader.filesize_ <= YON_FOOTER_LENGTH) {
		std::cerr << utility::timestamp("ERROR") << "File is corrupted!" << std::endl;
		return false;
	}

	// Seek to start of footer
	this->mImpl->basic_reader.stream_.seekg((uint64_t)this->mImpl->basic_reader.filesize_ - YON_FOOTER_LENGTH);
	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR") << "Failed to seek in file!" << std::endl;
		return false;
	}
	this->mImpl->basic_reader.stream_ >> this->global_footer;

	// Validate the global footer.
	if (this->global_footer.Validate() == false) {
		std::cerr << utility::timestamp("ERROR") << "Failed to validate footer!" << std::endl;
		return false;
	}

	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR") << "Failed to read file!" << std::endl;
		return false;
	}

	// Seek to start of file
	this->mImpl->basic_reader.stream_.seekg(0);
	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR") << "Failed to rewind file!" << std::endl;
		return false;
	}

	// Load header
	//this->stream >> this->global_header;
	char magic_string[TACHYON_MAGIC_HEADER_LENGTH];
	this->mImpl->basic_reader.stream_.read(&magic_string[0], TACHYON_MAGIC_HEADER_LENGTH);
	if (strncmp(&magic_string[0], &TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH) != 0) {
		std::cerr << utility::timestamp("ERROR") << "Failed to validate Tachyon magic string!" << std::endl;
		return false;
	}

	uint32_t l_data   = 0;
	uint32_t l_c_data = 0;
	utility::DeserializePrimitive(l_data, this->mImpl->basic_reader.stream_);
	utility::DeserializePrimitive(l_c_data, this->mImpl->basic_reader.stream_);

	yon_buffer_t header_uncompressed(l_data + 1024);
	yon_buffer_t header_compressed(l_c_data + 1024); header_compressed.n_chars_ = l_c_data;

	this->mImpl->basic_reader.stream_.read(header_compressed.data(), l_c_data);

	if (!this->mImpl->codec_manager.zstd_codec.Decompress(header_compressed, header_uncompressed)) {
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	assert(header_uncompressed.size() == l_data);
	header_uncompressed >> this->global_header; // parse header from buffer

	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR") << "Failed to get header!" << std::endl;
		return false;
	}

	// Keep track of start of data offset in the byte stream.
	// We use this information if spawning copies of this
	// reader object.
	this->b_data_start = this->mImpl->basic_reader.stream_.tellg();

	// Keep track of start position
	const uint64_t return_pos = this->mImpl->basic_reader.stream_.tellg();
	this->mImpl->basic_reader.stream_.seekg(this->global_footer.offset_end_of_data);
	this->mImpl->basic_reader.stream_ >> this->index;
	this->mImpl->basic_reader.stream_ >> this->mImpl->checksums;
	this->mImpl->basic_reader.stream_.seekg(return_pos);

	return(this->mImpl->basic_reader.stream_.good());
}

bool VariantReader::open(const std::string& filename) {
	this->mImpl->basic_reader.filename_ = filename;
	this->settings.input = filename;
	if (settings.keychain_file.size()) {
		if (this->LoadKeychainFile() == false)
			return false;
	}
	return(this->open());
}

bool VariantReader::NextBlock() {
	if (this->CheckNextValid() == false) return false;

	// Reset and re-use
	this->variant_container.clear();

	// Read packed header and footer bytes from the byte-stream.
	if (!this->variant_container.ReadHeaderFooter(this->mImpl->basic_reader.stream_))
		return false;

	// Decompress the footer.
	if (!this->mImpl->codec_manager.zstd_codec.Decompress(this->variant_container.footer_support)) {
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
		return false;
	}
	// Inject the decompressed footer into the container.
	this->variant_container.footer_support.data_uncompressed >> this->variant_container.footer;

	// Attempts to read a YON block with the settings provided.
	if (!this->variant_container.read(this->mImpl->basic_reader.stream_,
	                                 this->block_settings,
	                                 this->global_header))
	{
		return false;
	}

	// encryption manager ascertainment
	if (this->variant_container.header.controller.any_encrypted) {
		if (this->keychain.size() == 0) {
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption_manager_type encryption_manager;
		if (!encryption_manager.Decrypt(this->variant_container, this->keychain)) {
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if (!this->mImpl->codec_manager.Decompress(this->variant_container)) {
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// Checks need to made to ascertain that any data actually exists
	// and that if data was encrypted that any data was successfully
	// decrypted.

	// All passed
	return true;
}

bool VariantReader::NextBlockRaw(void) {
	if (this->CheckNextValid() == false) return false;

	// Reset and re-use
	this->variant_container.clear();

	if (!this->variant_container.ReadHeaderFooter(this->mImpl->basic_reader.stream_))
		return false;

	if (!this->mImpl->codec_manager.zstd_codec.Decompress(this->variant_container.footer_support)) {
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
		return false;
	}
	this->variant_container.footer_support.data_uncompressed >> this->variant_container.footer;

	// Attempts to read a YON block with the settings provided
	if (!this->variant_container.read(this->mImpl->basic_reader.stream_, this->block_settings, this->global_header))
		return false;

	// All passed
	return true;
}

bool VariantReader::CheckNextValid(void) {
	// If the stream is faulty then return
	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if ((uint64_t)this->mImpl->basic_reader.stream_.tellg() == this->global_footer.offset_end_of_data)
		return false;

	return true;
}

yon1_vb_t VariantReader::ReturnBlock(void) {
	block_entry_type c;

	if (this->CheckNextValid() == false)
		return c;

	if (!c.ReadHeaderFooter(this->mImpl->basic_reader.stream_))
		return(c);

	if (!this->mImpl->codec_manager.zstd_codec.Decompress(c.footer_support)) {
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	c.footer_support.data_uncompressed >> c.footer;

	// Attempts to read a YON block with the settings provided
	if (!c.read(this->mImpl->basic_reader.stream_, this->block_settings, this->global_header))
		return(c);

	// encryption manager ascertainment
	if (c.header.controller.any_encrypted) {
		if (this->keychain.size() == 0) {
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return(c);
		}

		encryption_manager_type encryption_manager;
		if (!encryption_manager.Decrypt(c, this->keychain)) {
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return(c);
		}
	}

	// Internally decompress available data
	if (!this->mImpl->codec_manager.Decompress(c)) {
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return(c);
	}

	return(c);
}

bool VariantReader::GetBlock(const index_entry_type& index_entry) {
	// If the stream is not good then return.
	if (!this->mImpl->basic_reader.stream_.good()) {
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// Seek to target block id with the help of the linear index.
	this->mImpl->basic_reader.stream_.seekg(index_entry.byte_offset);
	if (this->mImpl->basic_reader.stream_.good() == false) {
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to seek to given offset using target index entry!" << std::endl;
		return(false);
	}

	// Load the next block if possible.
	return(this->NextBlock());
}

bool VariantReader::SeekBlock(const uint32_t& block_id) {
	const uint64_t offset = this->GetIndex()[block_id].byte_offset;
	this->mImpl->basic_reader.stream_.seekg(offset);
	if (this->mImpl->basic_reader.stream_.good() == false) {
		std::cerr << "failed to seek" << std::endl;
		return false;
	}
	return true;
}

bool VariantReader::LoadKeychainFile(void) {
	std::ifstream keychain_reader(settings.keychain_file, std::ios::binary | std::ios::in);
	if (!keychain_reader.good()) {
		std::cerr << tachyon::utility::timestamp("ERROR") <<  "Failed to open keychain: " << settings.keychain_file << "..." << std::endl;
		return false;
	}

	keychain_reader >> this->keychain;
	if (!keychain_reader.good()) {
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse keychain..." << std::endl;
		return false;
	}
	return true;
}

uint64_t VariantReader::OutputVcfLinear(void) {
	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	yon_buffer_t buf(100000);
	while (this->NextBlock()) {
		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		if (this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.gt_ppa);

		for (uint32_t i = 0; i < vc.size(); ++i) {
			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}
			vc[i].ToVcfString(this->global_header, buf, this->GetBlockSettings().display_static, this->gt_exp);
			std::cout.write(buf.data(), buf.size());
			buf.reset();
		}
	}

	return 0;
}

uint64_t VariantReader::OutputVcfSearch(void) {
	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;
	yon_buffer_t buf(100000);

	for (uint32_t i = 0; i < this->mImpl->interval_container.GetBlockList().size(); ++i) {
		this->GetBlock(this->mImpl->interval_container.GetBlockList()[i]);

		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		if (this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.gt_ppa);

		for (uint32_t i = 0; i < vc.size(); ++i) {
			if ((this->*filter_intervals)(vc[i]) == false)
				continue;

			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}

			vc[i].ToVcfString(this->global_header, buf, this->GetBlockSettings().display_static, this->gt_exp);
			std::cout.write(buf.data(), buf.size());
			buf.reset();
		}
	}

	return 0;
}

uint64_t VariantReader::OutputYonLinear(void) {
	VariantWriterInterface* writer = nullptr;
	if (settings.output == "-" || settings.output.size() == 0)
		writer = new VariantWriterStream();
	else
		writer = new VariantWriterFile();

	if (writer->open(settings.output) == false) {
		std::cerr << "failed to open: " << settings.output << std::endl;
		return false;
	}

	writer->n_s = this->GetHeader().GetNumberSamples();
	writer->index.Setup(this->GetHeader().contigs_);
	writer->settings = this->settings;
	std::cerr << "writer settigns is=" << writer->settings.checkpoint_n_snps << std::endl;
	writer->UpdateHeaderView(this->global_header, this->settings.GetSettingsString());
	writer->WriteFileHeader(this->global_header);

	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	while (this->NextBlock()) {
		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		if (this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.gt_ppa);

		for (uint32_t i = 0; i < vc.size(); ++i) {
			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}
			if (vc[i].gt != nullptr) vc[i].gt->Expand();
			*writer += vc[i];
		}
	}

	writer->close();
	delete writer;

	return 0;
}

uint64_t VariantReader::OutputYonSearch(void) {
	VariantWriterInterface* writer = nullptr;
	if (settings.output == "-" || settings.output.size() == 0)
		writer = new VariantWriterStream();
	else
		writer = new VariantWriterFile();

	if (writer->open(settings.output) == false) {
		std::cerr << "failed to open: " << settings.output << std::endl;
		return false;
	}

	writer->n_s = this->GetHeader().GetNumberSamples();
	writer->index.Setup(this->GetHeader().contigs_);
	writer->settings = this->settings;
	writer->UpdateHeaderView(this->global_header, this->settings.GetSettingsString());
	writer->WriteFileHeader(this->global_header);

	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;

	for (uint32_t i = 0; i < this->mImpl->interval_container.GetBlockList().size(); ++i) {
		this->GetBlock(this->mImpl->interval_container.GetBlockList()[i]);

		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		if (this->GetBlockSettings().annotate_extra && this->settings.group_file.size())
			this->occ_table.BuildTable(this->variant_container.gt_ppa);

		for (uint32_t i = 0; i < vc.size(); ++i) {
			if ((this->*filter_intervals)(vc[i]) == false)
				continue;

			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}

			if (vc[i].gt != nullptr) vc[i].gt->Expand();
			*writer += vc[i];
		}
	}

	writer->close();
	delete writer;

	return 0;
}

uint64_t VariantReader::OutputRecords(void) {
	this->UpdateHeaderView();
	this->mImpl->interval_container.Build(this->global_header);

	// Load encryption keychain if available.
	if (this->settings.keychain_file.size())
		this->LoadKeychainFile();

	if (this->settings.annotate_genotypes)
		this->GetHeader().AddGenotypeAnnotationFields();

	// Occ
	if (this->settings.group_file.size()) {
		if (this->occ_table.ReadTable(this->settings.group_file, this->GetHeader(), '\t') == false) {
			return(0);
		}

		this->global_header.AddGenotypeAnnotationFields(this->occ_table.row_names);
	}

	// Any htslib
	if (this->settings.use_htslib) {
		if (this->mImpl->interval_container.size()) return(this->OutputHtslibVcfSearch());
		else return(this->OutputHtslibVcfLinear());
	}

	// Tachyon
	if (this->settings.output_type == 'y') {
		if (this->mImpl->interval_container.size()) return(this->OutputYonSearch());
		else return(this->OutputYonLinear());
	}

	// Vcf
	if (this->GetBlockSettings().show_vcf_header)
		this->global_header.PrintVcfHeader(std::cout);

	if (this->mImpl->interval_container.size()) return(this->OutputVcfSearch());
	else return(this->OutputVcfLinear());
}

uint64_t VariantReader::OutputHtslibVcfLinear(void) {
	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	// Open a htslib file handle for the target output
	// destination.
	char hts_stream_type[2];
	hts_stream_type[0] = 'w';
	hts_stream_type[1] = this->settings.output_type;
	htsFile* fp = hts_open(this->settings.output.c_str(), hts_stream_type);

	// Add extra htslib compression threads if desired.
	int n_extra_threads = std::thread::hardware_concurrency();
	if (n_extra_threads) {
		int ret = hts_set_threads(fp, n_extra_threads);
		if (ret < 0) {
			std::cerr << "failed to open multiple handles" << std::endl;
			return 0;
		}
	}

	// Convert the internal yon header to a bcf_hdr_t
	// structure.
	bcf_hdr_t* hdr = this->mImpl->ConvertVcfHeader(this->GetHeader(), !this->settings.drop_format);
	if ( bcf_hdr_write(fp, hdr) != 0 ) {
		std::cerr << "Failed to write header to " << this->settings.output << std::endl;
		exit(1);
	}

	// Initialize an empty record that we will keep
	// reusing as we iterate over available yon records.
	bcf1_t *rec = bcf_init1();

	// Iterate over available blocks.
	while (this->NextBlock()) {
		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		if (this->GetBlockSettings().annotate_extra)
			this->occ_table.BuildTable(this->variant_container.gt_ppa);

		// Iterate over available records in this block.
		for (uint32_t i = 0; i < vc.size(); ++i) {
			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateOcc(this->occ_table);
				vc[i].EvaluateOccSummary(true);

				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
				vc[i].AddGenotypeStatisticsOcc(this->global_header, this->occ_table.row_names);
			}

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			vc[i].OutputHtslibVcfInfo(rec, hdr);
			vc[i].OutputHtslibVcfFormat(rec, hdr, this->GetBlockSettings().display_static & YON_BLK_BV_GT, this->gt_exp);
			vc[i].OutputHtslibVcfFilter(rec, hdr);


			if ( bcf_write1(fp, hdr, rec) != 0 ) {
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

uint64_t VariantReader::OutputHtslibVcfSearch(void) {
	if (this->gt_exp == nullptr)
		this->gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];

	// Open a htslib file handle for the target output
	// destination.
	char hts_stream_type[2];
	hts_stream_type[0] = 'w'; hts_stream_type[1] = this->settings.output_type;
	htsFile *fp = hts_open(this->settings.output.c_str(), hts_stream_type);

	// Add extra htslib compression threads if desired.
	int n_extra_threads = std::thread::hardware_concurrency();
	if (n_extra_threads) {
		int ret = hts_set_threads(fp, n_extra_threads);
		if (ret < 0) {
			std::cerr << "failed to open multiple handles" << std::endl;
			return 0;
		}
	}

	// Convert the internal yon header to a bcf_hdr_t
	// structure.
	bcf_hdr_t* hdr = this->mImpl->ConvertVcfHeader(this->GetHeader(), !this->settings.drop_format);
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
	while (this->NextBlock()) {
		yon1_vc_t vc(this->GetCurrentContainer(), this->global_header);

		// Iterate over available records in this block.
		for (uint32_t i = 0; i < vc.size(); ++i) {
			if ((this->*filter_intervals)(vc[i]) == false)
				continue;

			if (this->variant_filters.Filter(vc[i], i) == false)
				continue;

			if (this->GetBlockSettings().annotate_extra) {
				vc[i].EvaluateSummary(true);
				vc[i].AddGenotypeStatistics(this->global_header);
			}

			vc[i].UpdateHtslibVcfRecord(rec, hdr);
			vc[i].OutputHtslibVcfInfo(rec, hdr);
			vc[i].OutputHtslibVcfFormat(rec, hdr, this->GetBlockSettings().display_static & YON_BLK_BV_GT, this->gt_exp);
			vc[i].OutputHtslibVcfFilter(rec, hdr);

			if ( bcf_write1(fp, hdr, rec) != 0 ) {
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

bool VariantReader::FilterIntervals(const yon1_vnt_t& entry) const { return(this->mImpl->interval_container.FindOverlaps(entry).size()); }

bool VariantReader::AddIntervals(std::vector<std::string>& interval_strings) {
	return(this->mImpl->interval_container.ParseIntervals(interval_strings, this->global_header, this->index));
}

void VariantReader::UpdateHeaderView(void) {
	VcfExtra e;
	e.key = "tachyon_viewVersion";
	e.value = TACHYON_PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  TACHYON_PROGRAM_NAME + '-' + TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	this->GetHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetHeader().extra_fields_.push_back(e);
	e.key = "tachyon_viewCommand";
	e.value = LITERAL_COMMAND_LINE;
	this->GetHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetHeader().extra_fields_.push_back(e);
	e.key = "tachyon_viewCommandSettings";
	e.value = this->GetSettings().GetSettingsString();
	this->GetHeader().literals_ += "##" + e.key + "=" + e.value + '\n';
	this->GetHeader().extra_fields_.push_back(e);
}

/*
bool VariantReader::Benchmark(const uint32_t threads) {
	std::cerr << utility::timestamp("LOG") << "Starting benchmark with " << threads << " threads..." << std::endl;
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::LoadData);
	this->mImpl->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::UncompressData);
	this->mImpl->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::EvaluateData);
	this->mImpl->basic_reader.stream_.close();
	this->BenchmarkWrapper(threads, &VariantSlavePerformance::EvaluateRecords);
	this->mImpl->basic_reader.stream_.close();
	return(true);
}

bool VariantReader::BenchmarkWrapper(const uint32_t threads, bool(VariantSlavePerformance::*func)(yon1_vb_t*&)) {
	if (!this->open(settings.input)) {
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	this->GetBlockSettings().LoadAll(true);

	// Begin
	algorithm::Timer timer; timer.Start();
	const uint32_t n_threads = threads;
	yon_producer_vblock<VariantReader> prd(n_threads);
	prd.Setup(&VariantReader::NextBlockRaw, *this, this->variant_container);
	prd.Start();

	yon_consumer_vblock<VariantSlavePerformance>* csm = new yon_consumer_vblock<VariantSlavePerformance>[n_threads];
	VariantSlavePerformance* slave = new VariantSlavePerformance[n_threads];

	for (int i = 0; i < n_threads; ++i) {
		csm[i].thread_id      = i;
		csm[i].data_available = &prd.data_available;
		csm[i].data_pool      = &prd.data_pool;

		slave[i].settings = this->GetBlockSettings();
		slave[i].gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];
		slave[i].global_header = &this->global_header;

		csm[i].Start(func, slave[i]);
	}
	// Join consumer and producer threads.
	for (uint32_t i = 0; i < n_threads; ++i) csm[i].thread_.join();

	prd.all_finished = true;
	prd.thread_.join();

	for (uint32_t i = 1; i < n_threads; ++i) slave[0] += slave[i];

	delete [] slave;
	delete [] csm;

	const double time_elapsed = timer.Elapsed().count();
	std::cerr << utility::timestamp("LOG") << time_elapsed << "\t" << slave[0].data_loaded << "\t" << (double)slave[0].data_loaded/time_elapsed/1e6 << "\t" << slave[0].data_uncompressed << "\t" << (double)slave[0].data_uncompressed/time_elapsed/1e6 << std::endl;
	return true;
}
*/

bool VariantReader::Stats(void) {
	const uint32_t n_threads = std::thread::hardware_concurrency();
	//const uint32_t n_threads = 1;
	yon_producer_vblock<VariantReader> prd(n_threads);
	prd.Setup(&VariantReader::NextBlockRaw, *this, this->variant_container);
	prd.Start();

	yon_consumer_vblock<VariantSlaveTsTv>* csm = new yon_consumer_vblock<VariantSlaveTsTv>[n_threads];

	VariantSlaveTsTv* slave_test = new VariantSlaveTsTv[n_threads];

 	for (int i = 0; i < n_threads; ++i) {
		csm[i].thread_id      = i;
		csm[i].data_available = &prd.data_available;
		csm[i].data_pool      = &prd.data_pool;

		slave_test[i].settings = this->GetBlockSettings();
		slave_test[i].s.SetSize(this->global_header.GetNumberSamples());
		slave_test[i].s_local.SetSize(this->global_header.GetNumberSamples());
		slave_test[i].gt_exp = new yon_gt_rcd[this->global_header.GetNumberSamples()];
		slave_test[i].global_header = &this->global_header;

		csm[i].Start(&VariantSlaveTsTv::GatherGenotypeStatistics, slave_test[i]);
	}
	// Join consumer and producer threads.
	for (uint32_t i = 0; i < n_threads; ++i) csm[i].thread_.join();

	prd.all_finished = true;
	prd.thread_.join();

	yon_stats_tstv s(this->GetHeader().GetNumberSamples());
	for (uint32_t i = 0; i < n_threads; ++i) s += slave_test[i].s;

	delete [] slave_test;
	delete [] csm;

	yon_buffer_t json_buffer(250000);
	s.Evaluate();
	s.ToJsonString(json_buffer, this->global_header.samples_);

	std::cout.write(json_buffer.data(), json_buffer.size());
	std::cout.flush();

	return true;
}

bcf_hdr_t* VariantReader::VariantReaderImpl::ConvertVcfHeaderLiterals(const yon_vnt_hdr_t& hdr, const bool add_format) {
	std::string internal = hdr.literals_;
	internal += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (hdr.samples_.size() && add_format) {
		internal += "\tFORMAT\t";
		internal += hdr.samples_[0];
		for (size_t i = 1; i < hdr.samples_.size(); ++i)
			internal += "\t" + hdr.samples_[i];
	}
	internal += "\n";

	bcf_hdr_t* bcf_hdr = bcf_hdr_init("r");
	int ret = bcf_hdr_parse(bcf_hdr, (char*)internal.c_str());
	if (ret != 0) {
		std::cerr << "failed to get bcf header from literals" << std::endl;
		bcf_hdr_destroy(bcf_hdr);
		return(nullptr);
	}

	return(bcf_hdr);
}

bcf_hdr_t* VariantReader::VariantReaderImpl::ConvertVcfHeader(const yon_vnt_hdr_t& hdr, const bool add_format) {
	std::string internal = hdr.ToString(true);
	internal += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (hdr.samples_.size() && add_format) {
		internal += "\tFORMAT\t";
		internal += hdr.samples_[0];
		for (size_t i = 1; i < hdr.samples_.size(); ++i)
			internal += "\t" + hdr.samples_[i];
	}
	internal += "\n";

	bcf_hdr_t* bcf_hdr = bcf_hdr_init("r");
	int ret = bcf_hdr_parse(bcf_hdr, (char*)internal.c_str());
	if (ret != 0) {
		std::cerr << "failed to get bcf header from literals" << std::endl;
		bcf_hdr_destroy(bcf_hdr);
		return(nullptr);
	}

	return(bcf_hdr);
}

}
