#ifndef CORE_VARIANT_IMPORT_SETTINGS_H_
#define CORE_VARIANT_IMPORT_SETTINGS_H_

#include <thread>

namespace tachyon{

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

}

#endif /* CORE_VARIANT_IMPORT_SETTINGS_H_ */
