/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_VARIANT_IMPORTER_H_
#define TACHYON_VARIANT_IMPORTER_H_

#include <string>
#include <cstdint>
#include <memory>

namespace tachyon {

struct VariantImporterSettings {
public:
	VariantImporterSettings();
	~VariantImporterSettings() = default;

	std::string GetInterpretedString(void) const;

	inline void SetVerbose(const bool yes = true) { this->verbose = yes; }
	inline void SetInputFile(const std::string& input_name) { this->input_file = input_name; }
	inline void SetOutputPrefix(const std::string& output_prefix) { this->output_prefix = output_prefix; }
	inline void SetThreads(const int32_t n_threads) { this->n_threads = n_threads; }
	inline void SetHtslibThreads(const int32_t n_threads) { this->htslib_extra_threads = n_threads; }
	inline void SetPermute(const bool yes) { this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes) { this->encrypt_data = yes; }
	inline void SetCompressionLevel(const int32_t compression_level) { this->compression_level = compression_level; }
	inline void SetCheckpointBases(const int32_t bases) { this->checkpoint_bases = bases; }
	inline void SetCheckpointVariants(const int32_t variants) { this->checkpoint_n_snps = variants; }
	inline void SetPermuteGenotypes(const bool yes = true) { this->permute_genotypes = yes; }

public:
	bool verbose;
	bool permute_genotypes; // permute GT flag
	bool encrypt_data; // encryption flag
	int32_t checkpoint_n_snps; // number of variants until checkpointing
	int32_t checkpoint_bases; // number of bases until checkpointing
	int32_t n_threads; // number of parallel importer threads
	int32_t compression_level; // compression level sent to ZSTD
	int32_t htslib_extra_threads; // extra threads for compress/decompress htslib
	int32_t info_end_key; // key mapping to the INFO field END
	int32_t info_svlen_key; // key mapping to the INFO field SVLEN
	std::string input_file; // input file name
	std::string output_prefix; // output file prefix
};

class VariantImporter {
public:
	typedef VariantImporter           self_type;
	typedef VariantImporterSettings   settings_type;

public:
	VariantImporter();
	VariantImporter(const settings_type& settings);
	VariantImporter(const self_type& other) = delete;
	VariantImporter& operator=(const self_type& other) = delete;
	VariantImporter(self_type&& other) = delete;
	VariantImporter& operator=(self_type&& other) = delete;
	~VariantImporter();

	/**<
	 * Imports a Vcf file into a Tachyon archive using the user-provided settings
	 * struct.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build();

	inline settings_type& GetSettings(void) { return(this->settings_); }
	inline const settings_type& GetSettings(void) const { return(this->settings_); }

private:
	settings_type settings_;

	// Pimpl idiom
	class VariantImporterImpl;
	std::unique_ptr<VariantImporterImpl> mImpl;
};


} /* namespace Tachyon */

#endif /* VARIANT_IMPORTER_H_ */
