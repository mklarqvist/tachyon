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
#ifndef TACHYON_VARIANT_READER_H_
#define TACHYON_VARIANT_READER_H_

#include <regex>
#include <cmath>
#include <thread>

#include "encryption.h"
#include "primitive_container.h"
#include "header_footer.h"
#include "variant_reader_filters.h"
#include "index.h"
#include "support_vcf.h"
#include "variant_container.h"

namespace tachyon {

struct VariantReaderSettings {
public:
	typedef VariantReaderSettings self_type;

public:
	VariantReaderSettings();
	~VariantReaderSettings() = default;

	/**<
	 * Construct a string with the internal interpreted parameters
	 * @return Returns a string
	 */
	std::string GetSettingsString(void) const;

	inline void SetThreads(const int32_t n_threads){ this->n_threads = n_threads; }
	inline void SetPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void SetCompressionLevel(const int32_t compression_level){ this->compression_level = compression_level; }
	inline void SetCheckpointBases(const int32_t bases){ this->checkpoint_bases = bases; }
	inline void SetCheckpointVariants(const int32_t variants){ this->checkpoint_n_snps = variants; }
	inline void SetPermuteGenotypes(const bool yes = true){ this->permute_genotypes = yes; }

public:
	bool drop_format, header_only, show_header, annotate_genotypes;
	bool use_htslib;
	std::string input, output, group_file, keychain_file;
	char output_type;

	// If writing to a yon archive
	bool permute_genotypes; // permute GT flag
	bool encrypt_data; // encryption flag
	int32_t checkpoint_n_snps; // number of variants until checkpointing
	int32_t checkpoint_bases; // number of bases until checkpointing
	int32_t n_threads; // number of parallel importer threads
	int32_t compression_level; // compression level sent to ZSTD
};

class VariantReader {
public:
	typedef VariantReader         self_type;
	typedef yon_buffer_t          buffer_type;
	typedef yon_vnt_hdr_t         header_type;
	typedef yon_ftr_t             footer_type;
	typedef yon_vb_settings       block_settings_type;
	typedef VariantReaderSettings settings_type;
	typedef yon_index_t           index_type;
	typedef yon1_idx_rec          index_entry_type;
	typedef yon1_vb_t             block_entry_type;
	typedef VariantReaderFilters  variant_filter_type;
	typedef EncryptionDecorator   encryption_manager_type;
	typedef Keychain              keychain_type;

private:
	// Function pointer to interval slicing.
	typedef bool (self_type::*filter_intervals_function)(const yon1_vnt_t& rcd) const;

public:
	VariantReader();
	VariantReader(const std::string& filename);
	virtual ~VariantReader();

	/**<
	 * Retrieve current settings records. Settings object is used prior to
	 * each read of a `block` and can be freely modified before each call.
	 * Modifying the settings after a read block has been invoked has no
	 * effect on the loaded `block` data.
	 * @return A reference instance of the block settings object
	 */
	inline block_settings_type& GetBlockSettings(void){ return(this->block_settings); }
	inline const block_settings_type& GetBlockSettings(void) const{ return(this->block_settings); }

	/**<
	 * Retrieve current settings for the variant reader. This settings
	 * object controls the parsing/output of the reader itself. This is
	 * unlike the `block_settings_type` that controls the `DataBlock`
	 * parsing.
	 * @return A reference instance of the settings object
	 */
	inline settings_type& GetSettings(void){ return(this->settings); }
	inline const settings_type& GetSettings(void) const{ return(this->settings); }

	/**<
	 * Retrieve the current filter settings for the variant reader. This
	 * object controls the pointers to filter applied to each variant.
	 * @return A reference instance of the filter object
	 */
	inline variant_filter_type& GetFilterSettings(void){ return(this->variant_filters); }

	// Basic accessors
	inline header_type& GetHeader(void){ return(this->global_header); }
	inline const header_type& GetHeader(void) const{ return(this->global_header); }
	inline footer_type& GetFooter(void){ return(this->global_footer); }
	inline const footer_type& GetFooter(void) const{ return(this->global_footer); }
	inline index_type& GetIndex(void){ return(this->index); }
	inline const index_type& GetIndex(void) const{ return(this->index); }
	inline block_entry_type& GetCurrentContainer(void){ return(this->variant_container); }
	inline const block_entry_type& GetCurrentContainer(void) const{ return(this->variant_container); }

	/**<
	 * Opens a YON file. Performs all prerequisite checks and loads all
	 * auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void);

	/**<
	 * Opens a YON file. Performs all prerequisite checks and loads all
	 * auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(const std::string& filename);

	/**<
	 * Overloaded operator for blocks. Useful when looping over a range of
	 * blocks. This happens frequently in parallel programs.
	 * @param index Block index value in range [0..n_blocks)
	 * @return      Returns TRUE if operation was successful or FALSE otherwise
	 */
	bool operator[](const uint32_t position);

	/**<
	 * Not implemented
	 * @param chromosome_name
	 * @return
	 */
	bool SeekToBlockChromosome(const std::string& chromosome_name);

	/**<
	 * Not implemented
	 * @param chromosome_name
	 * @param from_bp_position
	 * @param to_bp_position
	 * @return
	 */
	bool SeekToBlockChromosome(const std::string& chromosome_name,
	                           const uint32_t from_bp_position,
	                           const uint32_t to_bp_position);

	/**<
	 * Get the next YON block in-order. The NextBlockRaw() simply loads
	 * the appropriate data into memory without decrypting and uncompressing.
	 * The NextBlock() function performs these additional steps.
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool NextBlock(void);
	bool NextBlockRaw(void);

	/**<
	 * Predicate for checking if it is possible to retrieve another block.
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool CheckNextValid(void);

	block_entry_type ReturnBlock(void);

	/**<
	 * Get the target YON block that matches the provided index entry. Internally
	 * this function seeks to the offset described in the index entry and then
	 * invokes NextBlock().
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool GetBlock(const index_entry_type& index_entry);


	/**<
	 * Seeks to a specific YON block without loading anything. This allows the
	 * user to seek to a specific block and change the block_settings (i.e.
	 * what fields to load) and then invoke NextBlock() for example.
	 * @param block_id Target yon block id.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool SeekBlock(const uint32_t& block_id);

	/**<
	 * Open a yon encryption keychain and store the loaded object in this reader.
	 * @param path File path to archive.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadKeychainFile(void);

	/**<
	 * Wrapper function for iteratively constructing output data.
	 * Internally selects the appropriate subroutine to execute:
	 * 1) Using linear yon
	 * 2) Using non-linear yon
	 * 3) Using linear htslib
	 * 4) Using non-linear htslib
	 * @return Returns the number of variants processed.
	 */
	uint64_t OutputRecords(void);
	uint64_t OutputVcfLinear(void);
	uint64_t OutputVcfSearch(void);
	uint64_t OutputYonLinear(void);
	uint64_t OutputYonSearch(void);
	uint64_t OutputHtslibVcfLinear(void);
	uint64_t OutputHtslibVcfSearch(void);

	// Filter interval intersection and dummy version
	inline bool FilterIntervalsDummy(const yon1_vnt_t& entry) const{ return true; }
	bool FilterIntervals(const yon1_vnt_t& entry) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool AddIntervals(std::vector<std::string>& interval_strings);

	/**<
	 * Adds provenance tracking information to the yon header. This data
	 * corresponds to the version of tachyon, versions of the linked
	 * libraries, the literal command line provided to the tachyon view
	 * subroutine, and the internal interpreted setting used in JSON
	 * format.
	 */
	void UpdateHeaderView(void);

	bool Stats(void);

private:
	// Pimpl idiom
	class VariantReaderImpl;
	std::unique_ptr<VariantReaderImpl> mImpl;

private:
	uint64_t                b_data_start;
	block_entry_type        variant_container;
	block_settings_type     block_settings;
	settings_type           settings;
	variant_filter_type     variant_filters;
	header_type             global_header;
	footer_type             global_footer;
	index_type              index;
	keychain_type           keychain;
	yon_occ                 occ_table;

	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	yon_gt_rcd* gt_exp;
};

}

#endif
