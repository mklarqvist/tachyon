#ifndef VARIANT_IMPORTER_H_
#define VARIANT_IMPORTER_H_

#include <unordered_map>

#include "io/variant_import_writer.h"

namespace tachyon {

struct VariantImporterSettings {
public:
	VariantImporterSettings();
	~VariantImporterSettings() = default;

	std::string GetInterpretedString(void) const;

	inline void SetInputFile(const std::string& input_name){ this->input_file = input_name; }
	inline void SetOutputPrefix(const std::string& output_prefix){ this->output_prefix = output_prefix; }
	inline void SetThreads(const uint32_t n_threads){ this->n_threads = n_threads; }
	inline void SetPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void SetCompressionLevel(const uint32_t compression_level){ this->compression_level = compression_level; }

public:
	bool permute_genotypes; // permute GT flag
	bool encrypt_data; // encryption flag
	uint32_t checkpoint_n_snps; // number of variants until checkpointing
	uint32_t checkpoint_bases; // number of bases until checkpointing
	uint32_t n_threads; // number of parallel importer threads
	int32_t info_end_key; // key mapping to the INFO field END
	int32_t info_svlen_key; // key mapping to the INFO field SVLEN
	uint32_t compression_level; // compression level sent to ZSTD
	std::string input_file;   // input file name
	std::string output_prefix;// output file prefix
	uint32_t htslib_extra_threads; // extra threads for compress/decompress htslib
};

class VariantImporter {
public:
	typedef VariantImporter           self_type;
	typedef VariantWriterInterface    writer_interface_type;
	typedef VariantWriterFile         writer_file_type;
	typedef VariantWriterStream       writer_stream_type;
	typedef VariantImporterSettings   settings_type;

public:
	VariantImporter();
	VariantImporter(const settings_type& settings);
	VariantImporter(const self_type& other) = delete;
	VariantImporter& operator=(const self_type& other) = delete;
	~VariantImporter();

	bool Build();
	void SetWriterTypeFile(void);
	void SetWriterTypeStream(void);

private:
	settings_type settings_;
	writer_interface_type* writer; // writer

	// Pimpl Idiom
	class VariantImporterImpl;
	std::unique_ptr<VariantImporterImpl> mImpl;
};


} /* namespace Tachyon */

#endif /* VARIANT_IMPORTER_H_ */
