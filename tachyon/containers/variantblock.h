#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include "../algorithm/permutation/permutation_manager.h"
#include "../core/variant_importer_container_stats.h"
#include "../io/vcf/VCFHeader.h"
#include "components/datablock_header.h"
#include "datablock_settings.h"
#include "datacontainer.h"
#include "../core/meta_entry.h"

namespace tachyon{
namespace containers{

/**
 * Primary Tachyon block object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class VariantBlock{
	typedef VariantBlock                    self_type;
	typedef DataContainer                   container_type;
	typedef algorithm::PermutationManager   permutation_type;
	typedef DataBlockHeader                 block_header_type;
	typedef DataBlockFooter                 block_footer_type;
	typedef HashContainer                   hash_container_type;
	typedef HashVectorContainer             hash_vector_container_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef core::DataBlockSettings         settings_type;
	typedef support::VariantImporterContainerStats  import_stats_type;
	typedef DataContainerHeader             offset_type;
	typedef tachyon::core::MetaEntry        meta_entry_type;

public:
	VariantBlock();
	VariantBlock(const U32 n_info_fields, const U32 n_format_fields);
	~VariantBlock();

	/**< @brief Resize base container buffer streams
	 * Internal use only
	 * @param s Size in bytes
	 */
	void resize(const U32 s);

	/**< @brief Recycle structure without releasing memory
	 * Internal use only: Clears data by resetting
	 * pointers and values without releasing and
	 * reallocating the memory
	 */
	void clear(void);

	inline const U32& size(void) const{ return(this->header.n_variants); }

	inline U32 AddFieldINFO(const U32 fieldID){ return(this->info_fields.setGet(fieldID)); }
	inline U32 AddFieldFORMAT(const U32 fieldID){ return(this->format_fields.setGet(fieldID)); }
	inline U32 AddFieldFILTER(const U32 fieldID){ return(this->filter_fields.setGet(fieldID)); }

	inline const bool checkPatternsINFO(const U64& hash_pattern) const{
		U32 mapID = 0;
		return(this->info_patterns.getRaw(hash_pattern, mapID));
	}

	inline const bool checkPatternsFORMAT(const U64& hash_pattern) const{
		U32 mapID = 0;
		return(this->format_patterns.getRaw(hash_pattern, mapID));
	}

	inline const bool checkPatternsFILTER(const U64& hash_pattern) const{
		U32 mapID = 0;
		return(this->filter_patterns.getRaw(hash_pattern, mapID));
	}

	inline void addPatternINFO(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->info_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


			for(size_t i = 0; i < pattern.size(); ++i){
				std::cerr << pattern[i] << std::endl;
			}
			exit(1);
		}
	}

	inline void addPatternFORMAT(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->format_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


			for(size_t i = 0; i < pattern.size(); ++i){
				std::cerr << pattern[i] << std::endl;
			}
			exit(1);
		}
	}

	inline void addPatternFILTER(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->filter_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			std::cerr << this->format_patterns.size() << "," << this->format_fields.size() << "\t" << this->info_patterns.size() << "," << this->format_fields.size() << "\t" << this->filter_patterns.size() << "," << this->filter_fields.size() << std::endl;


			for(size_t i = 0; i < pattern.size(); ++i){
				std::cerr << pattern[i] << std::endl;
			}
			exit(1);
		}
	}

	/**<
	 * Finalize this block before writing to disk
	 * @param info_values
	 * @param info_patterns
	 * @param filter_values
	 * @param filter_patterns
	 * @param format_values
	 * @param format_patterns
	 * @return
	 */
	inline bool finalize(void){
		this->footer.n_info_streams   = this->info_fields.size();
		this->footer.n_filter_streams = this->filter_fields.size();
		this->footer.n_format_streams = this->format_fields.size();
		this->allocateDiskOffsets(this->info_fields.size(), this->format_fields.size(), this->filter_fields.size());
		this->updateBaseContainers();
		this->updateContainerSet(containers::DataBlockFooter::INDEX_INFO);
		this->updateContainerSet(containers::DataBlockFooter::INDEX_FORMAT);
		this->footer.constructBitVector(containers::DataBlockFooter::INDEX_INFO,   this->info_fields,   this->info_patterns);
		this->footer.constructBitVector(containers::DataBlockFooter::INDEX_FILTER, this->filter_fields, this->filter_patterns);
		this->footer.constructBitVector(containers::DataBlockFooter::INDEX_FORMAT, this->format_fields, this->format_patterns);
		return true;
	}

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading data from disk. Data
	 * read in this way is not checked for integrity here.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool read(std::ifstream& stream, settings_type& settings);

	/**<
	 * Standard way of writing out a YON block.
	 * @param stream       Target output stream
	 * @param stats_basic  Tracking for basic containers
	 * @param stats_info   Tracking for INFO containers
	 * @param stats_format Tracking for FORMAT containers
	 * @return             Returns TRUE upon success or FALSE otherwise
	 */
	bool write(std::ofstream& stream, import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format);

	// Add a meta entry
	/**<
	 * Overloaded operator for adding a variant entry
	 * @param meta_entry
	 * @return
	 */
	bool operator+=(meta_entry_type& meta_entry);
	inline bool operator<<(meta_entry_type& meta_entry){ return(*this += meta_entry); }

private:
	/**<
	 * Wrapper function for INFO, FORMAT, and FILTER disk
	 * offsets. Internally generates virtual disk offsets
	 * into the buffer stream where a given byte stream
	 * begins. This allows for random access searches to the
	 * data in question.
	 * @param n_info_fields   Number of INFO fields in this block
	 * @param n_format_fields Number of FORMAT fields in this block
	 * @param n_filter_fields Number of FILTER fields in this block
	 */
	inline void allocateDiskOffsets(const U32& n_info_fields, const U32& n_format_fields, const U32& n_filter_fields){
		this->footer.allocateDiskOffsets(n_info_fields, n_format_fields, n_filter_fields);
	}

	/**<
	 * Wrapper function (indirection) for invoking the updateContainer
	 * function for either all INFO or FORMAT containers
	 * @param target Enum target for groups of containers
	 */
	inline void updateContainerSet(DataBlockFooter::INDEX_BLOCK_TARGET target){
		// Determine target
		switch(target){
		case(DataBlockFooter::INDEX_BLOCK_TARGET::INDEX_INFO)   :
			return(this->updateContainer(this->info_containers, this->footer.n_info_streams));
			break;
		case(DataBlockFooter::INDEX_BLOCK_TARGET::INDEX_FORMAT) :
			return(this->updateContainer(this->format_containers, this->footer.n_format_streams));
			break;
		default: std::cerr << "unknown target type" << std::endl; exit(1);
		}
	}

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used primitive type) for strides and data; if possible
	 */
	void updateBaseContainers(void);

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used word-size) for strides and data; if possible
	 *
	 * @param container Data container
	 * @param length    Iterator length
	 */
	void updateContainer(container_type* container, const U32& length);

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used word-size) for strides and data; if possible
	 *
	 * @param container Data container
	 * @param reormat   Reformat boolean
	 */
	void updateContainer(container_type& container, bool reformat = true);

	void deltaEncode(container_type& container);

	const U64 __determineCompressedSize(void) const;

	inline void __updateHeader(offset_type& offset, const container_type& container){
		const U32 global_key = offset.data_header.global_key;
		offset = container.header;
		offset.data_header.global_key = global_key;
	}

	inline void __updateHeader(offset_type& offset, const container_type& container, const U32& virtual_offset){
		const U32 global_key = offset.data_header.global_key;
		offset = container.header;
		offset.data_header.global_key = global_key;
		offset.data_header.offset = virtual_offset;
	}

public:
	block_header_type header;
	block_footer_type footer;
	permutation_type  ppa_manager;
	container_type    meta_contig_container;
	container_type    meta_positions_container;
	container_type    meta_refalt_container;
	container_type    meta_controller_container;
	container_type    meta_quality_container;
	container_type    meta_names_container;
	container_type    meta_alleles_container;
	container_type    meta_info_map_ids;
	container_type    meta_format_map_ids;
	container_type    meta_filter_map_ids;
	container_type    gt_support_data_container; // data (1: diploid-rle, 2: diploid-other, 3: diploid-bcf, 4: other-ploidy-bcf), strides (n_objects OR ploidy for case 4)
	container_type    gt_rle8_container;
	container_type    gt_rle16_container;
	container_type    gt_rle32_container;
	container_type    gt_rle64_container;
	container_type    gt_simple_container;
	container_type*   info_containers;
	container_type*   format_containers;

	// Use during construction
	hash_container_type        info_fields;
	hash_container_type        format_fields;
	hash_container_type        filter_fields;
	hash_vector_container_type info_patterns;
	hash_vector_container_type format_patterns;
	hash_vector_container_type filter_patterns;

public:
	// Utility
	//size_t n_capacity_info_;
	//size_t n_capacity_format_;
	//U64    disk_offset_;       // utility primitive to support the construction of iterators over blocks
	container_type footer_support; // used internally only
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
