#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include "../algorithm/permutation/permutation_manager.h"
#include "../core/variant_importer_container_stats.h"
#include "../io/vcf/VCFHeader.h"
#include "components/datablock_header.h"
#include "datablock_settings.h"
#include "datacontainer.h"

namespace tachyon{
namespace containers{

/**
 * Primary Tachyon block object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class VariantBlock{
	typedef VariantBlock                   self_type;
	typedef DataContainer                  container_type;
	typedef algorithm::PermutationManager  permutation_type;
	typedef DataBlockHeader                header_type;
	typedef HashContainer                  hash_container_type;
	typedef HashVectorContainer            hash_vector_container_type;
	typedef io::BasicBuffer                buffer_type;
	typedef core::DataBlockSettings        settings_type;
	typedef support::VariantImporterContainerStats  import_stats_type;

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
		this->header.allocateDiskOffsets(n_info_fields, n_format_fields, n_filter_fields);
	}

	/**<
	 * Wrapper function (indirection) for invoking the updateContainer
	 * function for either all INFO or FORMAT containers
	 * @param target Enum target for groups of containers
	 */
	inline void updateContainerSet(DataBlockHeader::INDEX_BLOCK_TARGET target){
		// Determine target
		switch(target){
		case(DataBlockHeader::INDEX_BLOCK_TARGET::INDEX_INFO)   :
			return(this->updateContainer(this->info_containers, this->header.n_info_streams));
			break;
		case(DataBlockHeader::INDEX_BLOCK_TARGET::INDEX_FORMAT) :
			return(this->updateContainer(this->format_containers, this->header.n_format_streams));
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

	/**< @brief Updates the local (relative to block start) virtual file offsets
	 *  Internal use only (import): This functions updates
	 *  relative (virtual) file offsets into the byte stream
	 *  where a target container/object begins. This update
	 *  allows random access by having the block-header only
	 *
	 *  Each digital object must have an associated getDiskSize()
	 *  function that returns its byte length.
	 */
	void updateOffsets(void);

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

private:
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

public:
	header_type       header;
	permutation_type  ppa_manager;
	container_type    meta_chromosome_container;
	container_type    meta_position_container;
	container_type    meta_hot_container;
	container_type    meta_info_map_ids;
	container_type    meta_format_map_ids;
	container_type    meta_filter_map_ids;
	container_type    meta_cold_container;
	container_type    gt_support_data_container; // data (1: diploid-rle, 2: diploid-other, 3: diploid-bcf, 4: other-ploidy-bcf), strides (n_objects OR ploidy for case 4)
	container_type    gt_rle_container;
	container_type    gt_simple_container;
	container_type*   info_containers;
	container_type*   format_containers;

public:
	// Utility
	size_t n_capacity_info_;
	size_t n_capacity_format_;
	U64    disk_offset_;       // utility primitive to support the construction of iterators over blocks
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
