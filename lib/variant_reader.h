#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <regex>
#include <cmath>
#include <thread>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "algorithm/timer.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/interval_container.h"
#include "containers/meta_container.h"
#include "containers/primitive_group_container.h"
#include "containers/variant_block.h"
#include "containers/variant_block_container.h"
#include "core/footer/footer.h"
#include "core/header/variant_header.h"
#include "core/variant_reader_filters.h"
#include "core/variant_reader_objects.h"
#include "core/variant_reader_settings.h"
#include "core/data_block_settings.h"
#include "index/index.h"
#include "math/basic_vector_math.h"
#include "math/fisher_math.h"
#include "math/square_matrix.h"
#include "utility/support_vcf.h"
#include "io/basic_reader.h"

#include "core/ts_tv_object.h"

namespace tachyon{

class VariantReader{
public:
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef VariantHeader                          header_type;
	typedef core::Footer                           footer_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef DataBlockSettings                      block_settings_type;
	typedef VariantReaderSettings                  settings_type;
	typedef index::Index                           index_type;
	typedef index::VariantIndexEntry               index_entry_type;
	typedef algorithm::VariantDigestManager        checksum_type;
	typedef Keychain                               keychain_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;
	typedef containers::IntervalContainer          interval_container_type;
	typedef containers::VariantBlockContainer      variant_container_type;
	typedef VariantReaderFilters                   variant_filter_type;
	typedef algorithm::Interval<uint32_t, int64_t> interval_type;
	typedef io::BasicReader                        basic_reader_type;
	typedef EncryptionDecorator                    encryption_manager_type;

private:
	// Function pointer to interval slicing.
	typedef bool (self_type::*filter_intervals_function)(const meta_entry_type& meta_entry) const;

public:
	VariantReader();
	VariantReader(const std::string& filename);
	VariantReader(const self_type& other);
	virtual ~VariantReader();

	/**<
	 * Retrieve current settings records. Settings object is used
	 * prior to each read of a `block` and can be freely modified
	 * before each call. Modifying the settings after a read block
	 * has been invoked has no effect on the loaded `block` data.
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
	inline header_type& GetGlobalHeader(void){ return(this->global_header); }
	inline const header_type& GetGlobalHeader(void) const{ return(this->global_header); }
	inline footer_type& GetGlobalFooter(void){ return(this->global_footer); }
	inline const footer_type& GetGlobalFooter(void) const{ return(this->global_footer); }
	inline index_type& GetIndex(void){ return(this->index); }
	inline const index_type& GetIndex(void) const{ return(this->index); }
	inline size_t GetFilesize(void) const{ return(this->basic_reader.filesize_); }
	inline variant_container_type& GetCurrentContainer(void){ return(this->variant_container); }
	inline const variant_container_type& GetCurrentContainer(void) const{ return(this->variant_container); }

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void);

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	inline bool open(const std::string& filename){
		this->basic_reader.filename_ = filename;
		this->settings.input = filename;
		if(settings.keychain_file.size()){
			if(this->LoadKeychainFile() == false)
				return false;
		}
		return(this->open());
	}

	/**<
	 * Overloaded operator for blocks. Useful when
	 * looping over a range of blocks. This happens
	 * frequently in parallel programs.
	 * @param index Block index value in range [0..n_blocks)
	 * @return      Returns TRUE if operation was successful or FALSE otherwise
	 */
	bool operator[](const uint32_t position);

	/**<
	 * Not implemented
	 * @param position
	 * @return
	 */
	bool SeektoBlock(const uint32_t position);

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
	bool SeekToBlockChromosome(const std::string& chromosome_name, const uint32_t from_bp_position, const uint32_t to_bp_position);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool NextBlock(void);

	bool CheckNextValid(void){
		// If the stream is faulty then return
		if(!this->basic_reader.stream_.good()){
			std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
			return false;
		}

		// If the current position is the EOF then
		// exit the function
		if((uint64_t)this->basic_reader.stream_.tellg() == this->global_footer.offset_end_of_data)
			return false;

		return true;
	}

	containers::VariantBlockContainer ReturnBlock(void);

	/**<
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool GetBlock(const index_entry_type& index_entry);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the block_settings (i.e. what fields to load) and
	 * then invoke NextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool SeekBlock(const uint32_t& blockID){
		const uint64_t offset = this->GetIndex().GetLinearIndex().at(blockID).byte_offset;
		std::cerr << "offset is " << offset << std::endl;
		this->basic_reader.stream_.seekg(offset);
		if(this->basic_reader.stream_.good() == false){
			std::cerr << "failed to seek" << std::endl;
			return false;
		}
		return true;
	}

	/**<
	 * Open a yon encryption keychain and store the loaded object
	 * in this reader.
	 * @param path File path to archive.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool LoadKeychainFile(void);

	uint64_t OutputRecords(void);
	uint64_t OutputVcfLinear(void);
	uint64_t OutputVcfSearch(void);
	void OuputVcfWrapper(io::BasicBuffer& output_buffer, yon1_t& entry) const;
	void OutputInfoVcf(io::BasicBuffer& output_buffer, yon1_t& entry) const;
	void OutputFormatVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;
	void OutputFilterVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;

	uint64_t OutputHtslibVcfLinear(void);
	uint64_t OutputHtslibVcfSearch(void);
	void OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr, yon1_t& entry) const;
	void OutputHtslibVcfFormat(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_t& entry) const;
	void OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr, const yon1_t& entry) const;

	// Filter interval intersection and dummy version
	inline bool FilterIntervalsDummy(const meta_entry_type& meta_entry) const{ return true; }
	inline bool FilterIntervals(const meta_entry_type& meta_entry) const{ return(this->interval_container.FindOverlaps(meta_entry).size()); }

	// Calculations
	TACHYON_VARIANT_CLASSIFICATION_TYPE ClassifyVariant(const meta_entry_type& meta, const uint32_t& allele) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	inline bool AddIntervals(std::vector<std::string>& interval_strings){
		return(this->interval_container.ParseIntervals(interval_strings, this->global_header, this->index));
	}

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
	uint64_t b_data_start;
	basic_reader_type       basic_reader;
	variant_container_type  variant_container;
	block_settings_type     block_settings;
	settings_type           settings;
	variant_filter_type     variant_filters;
	header_type             global_header;
	footer_type             global_footer;
	index_type              index;
	checksum_type           checksums;
	codec_manager_type      codec_manager;
	keychain_type           keychain;
	interval_container_type interval_container;
	yon_occ                 occ_table;
};

}

#endif /* CORE_TACHYON_READER_H_ */
