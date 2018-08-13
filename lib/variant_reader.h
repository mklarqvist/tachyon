#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <regex>
#include <cmath>

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

namespace tachyon{

class VariantReader{
private:
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef VariantHeader                          header_type;
	typedef core::Footer                           footer_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef DataBlockSettings                      block_settings_type;
	typedef VariantReaderSettings                  settings_type;
	typedef index::Index                           index_type;
	typedef index::IndexEntry                      index_entry_type;
	typedef algorithm::VariantDigestManager        checksum_type;
	typedef encryption::Keychain<>                 keychain_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;
	typedef containers::IntervalContainer          interval_container_type;
	typedef containers::VariantBlockContainer      variant_container_type;
	typedef VariantReaderFilters                   variant_filter_type;
	typedef algorithm::Interval<U32, S64>          interval_type;
	typedef io::BasicReader                        basic_reader_type;
	typedef encryption::EncryptionDecorator        encryption_manager_type;

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
	bool operator[](const U32 position);

	/**<
	 * Not implemented
	 * @param position
	 * @return
	 */
	bool SeektoBlock(const U32 position);

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
	bool SeekToBlockChromosome(const std::string& chromosome_name, const U32 from_bp_position, const U32 to_bp_position);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool NextBlock(void);

	/**<
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool GetBlock(const index_entry_type& index_entry);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the block_settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool SeekBlock(const U32& blockID);

	U64 OutputVcf(void);
	U64 OutputVcfLinear(void);
	U64 OutputVcfSearch(void);
	void OuputVcfWrapper(io::BasicBuffer& output_buffer, yon1_t& entry) const;
	void OutputInfoVcf(io::BasicBuffer& output_buffer, yon1_t& entry) const;
	void OutputFormatVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;
	void OutputFilterVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;

	/**<
	 * Wrapper function to call internal functions `outputCustom` or `outputBlockVCF`.
	 * Decides internally what function to invoke.
	 * @return
	 */
	U64 outputVCF(void);

	// Filter interval intersection and dummy version
	inline bool FilterIntervalsDummy(const meta_entry_type& meta_entry) const{ return true; }
	inline bool FilterIntervals(const meta_entry_type& meta_entry) const{ return(this->interval_container.find_overlaps(meta_entry).size()); }

	// Calculations
	TACHYON_VARIANT_CLASSIFICATION_TYPE ClassifyVariant(const meta_entry_type& meta, const U32& allele) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	inline bool AddIntervals(std::vector<std::string>& interval_strings){
		return(this->interval_container.parseIntervals(interval_strings, this->global_header, this->index));
	}


	/**<
	 *
	 * @param path
	 * @return
	 */
	bool LoadKeychainFile(void){
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

	/**<
	 *
	 * @param stream
	 */
	void PrintHeaderVCF(std::ostream& stream = std::cout){
		this->global_header.literals_ += "##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		this->global_header.literals_ += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
		                             +   SSLeay_version(SSLEAY_VERSION) + ","
		                             +  "ZSTD-" + ZSTD_versionString()
		                             +  "; timestamp=" + tachyon::utility::datetime() + "\n";

		this->global_header.literals_ += "##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + "\n";
		this->global_header.literals_ += this->GetSettings().get_settings_string();
		this->global_header.literals_ += '\n';

		this->global_header.PrintVcfHeader(stream);
	}


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	/*
			if(target_flag_set & 512){
				// Classify
				buffer += ";VT=";
				buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->ClassifyVariant(objects.meta_container->at(position), 1)];

				for(U32 p = 2; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ',';
					buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->ClassifyVariant(objects.meta_container->at(position), p)];
				}
			}
			*/



private:
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
