#ifndef VARIANT_IMPORTER_H_
#define VARIANT_IMPORTER_H_

#include <unordered_map>

#include "algorithm/compression/compression_manager.h"
#include "algorithm/compression/genotype_encoder.h"
#include "algorithm/timer.h"
#include "containers/variant_block.h"
#include "core/variant_import_writer.h"
#include "core/variant_importer_container_stats.h"
#include "index/index_entry.h"
#include "index/index_index_entry.h"
#include "io/vcf_utils.h"
#include "support/helpers.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "core/footer/footer.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "algorithm/permutation/genotype_sorter.h"

namespace tachyon {

/**<
 * Settings for VariantImporter class
 */
struct VariantImporterSettings{
public:
	VariantImporterSettings() :
		permute_genotypes(true),
		encrypt_data(false),
		drop_invariant_sites(false),
		checkpoint_n_snps(1000),
		checkpoint_bases(10e6),
		n_threads(std::thread::hardware_concurrency()),
		info_end_key(-1),
		info_svlen_key(-1),
		compression_level(6),
		htslib_extra_threads(0)
	{}

	~VariantImporterSettings() = default;

	std::string GetInterpretedString(void) const{
		return(std::string("{\"input_file\":\"" + this->input_file +
		   "\",\"output_prefix\":\"" + this->output_prefix +
		   "\",\"checkpoint_snps\":" + std::to_string(this->checkpoint_n_snps) +
		   ",\"checkpoint_bases\":" + std::to_string(this->checkpoint_bases) +
		   ",\"compression_level\":" + std::to_string(this->compression_level) +
		   "}"
		));
	}

	inline void SetInputFile(const std::string& input_name){ this->input_file = input_name; }
	inline void SetOutputPrefix(const std::string& output_prefix){ this->output_prefix = output_prefix; }
	inline void SetThreads(const uint32_t n_threads){ this->n_threads = n_threads; }
	inline void SetPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void SetCompressionLevel(const uint32_t compression_level){ this->compression_level = compression_level; }

public:
	bool permute_genotypes;   // permute GT flag
	bool encrypt_data;        // encryption flag
	bool drop_invariant_sites;// drop sites that are invariant
	uint32_t checkpoint_n_snps;    // number of variants until checkpointing
	uint32_t checkpoint_bases;     // number of bases until checkpointing
	uint32_t n_threads;            // number of parallel importer threads
	int32_t info_end_key;         // key mapping to the INFO field END
	int32_t info_svlen_key;       // key mapping to the INFO field SVLEN
	uint32_t compression_level;    // compression level sent to ZSTD
	std::string input_file;   // input file name
	std::string output_prefix;// output file prefix
	uint32_t htslib_extra_threads;
};

class VariantImporter {
private:
	typedef VariantImporter                 self_type;
	typedef VariantImportWriterInterface    writer_interface_type;
	typedef VariantImportWriterFile         writer_file_type;
	typedef VariantImportWriterStream       writer_stream_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef index::IndexEntry               index_entry_type;
	typedef io::VcfReader                   vcf_reader_type;
	typedef containers::VcfContainer        vcf_container_type;
	typedef algorithm::CompressionManager   compression_manager_type;
	typedef algorithm::GenotypeSorter       radix_sorter_type;
	typedef algorithm::GenotypeEncoder      gt_encoder_type;
	typedef containers::DataContainer       stream_container;
	typedef containers::VariantBlock        block_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef core::MetaEntry                 meta_type;
	typedef VariantImporterSettings         settings_type;
	typedef std::unordered_map<uint32_t, uint32_t>    reorder_map_type;
	typedef std::unordered_map<uint64_t, uint32_t>    hash_map_type;

public:
	VariantImporter();
	VariantImporter(const settings_type& settings);
	~VariantImporter();

	bool Build();
	bool BuildParallel();

	void clear(void);

	inline void SetWriterTypeFile(void){ this->writer = new writer_file_type; }
	inline void SetWriterTypeStream(void){ this->writer = new writer_stream_type; }

private:
	bool BuildVCF();
	bool AddRecords(const vcf_container_type& container);
	bool AddRecord(const vcf_container_type& container, const uint32_t position, meta_type& meta);
	bool AddVcfInfo(const bcf1_t* record, meta_type& meta);
	bool AddVcfFormatInfo(const bcf1_t* record, meta_type& meta);
	bool AddVcfFilterInfo(const bcf1_t* record, meta_type& meta);
	bool IndexRecord(const bcf1_t* record, const meta_type& meta);
	bool AddVcfInfoPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddVcfFormatPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddVcfFilterPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddGenotypes(const vcf_container_type& container, meta_type* meta_entries);
	bool UpdateIndex();
	bool WriteBlock();
	bool WriteFinal(algorithm::VariantDigestManager& checksums);
	bool WriteKeychain(const encryption::Keychain<>& keychain);
	bool WriteYonHeader();
	void UpdateHeaderImport(VariantHeader& header);
	bool GenerateIdentifiers(void);

private:
	settings_type settings_; // internal settings
	bool GT_available_;

	// Stats
	import_stats_type stats_basic;
	import_stats_type stats_info;
	import_stats_type stats_format;

	// Read/write fields
	writer_interface_type* writer; // writer
	index_entry_type  index_entry; // streaming index entry
	radix_sorter_type permutator;  // GT permuter
	gt_encoder_type   encoder;     // RLE packer

	compression_manager_type compression_manager;

	// Data container
	block_type block;

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
	vcf_container_type vcf_container_;

	hash_map_type block_hash_map;
};


} /* namespace Tachyon */

#endif /* VARIANT_IMPORTER_H_ */
