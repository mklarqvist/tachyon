#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include "../algorithm/permutation/PermutationManager.h"
#include "../index/BlockIndex.h"
#include "Container.h"
#include "../io/vcf/VCFHeader.h"
#include "BlockEntrySettings.h"
#include "../core/ImporterStats.h"

namespace Tachyon{
namespace Core{

/**
 * Primary Tachyon object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class Block{
	typedef Block self_type;
	typedef Core::Container container_type;
	typedef Core::PermutationManager permutation_type;
	typedef Index::BlockIndex index_entry_type;
	typedef Core::Support::HashContainer hash_container_type;
	typedef Core::Support::HashVectorContainer hash_vector_container_type;
	typedef Index::BlockIndexOffsets offset_type;
	typedef Index::BlockIndexOffsetsHeader offset_minimal_type;
	typedef IO::BasicBuffer buffer_type;
	typedef BlockEntrySettings settings_type;
	typedef Tachyon::Support::ImporterStats import_stats_type;

public:
	Block();
	~Block();

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

	inline const U32 size(void) const{ return(this->index_entry.n_variants); }

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
		this->index_entry.allocateDiskOffsets(n_info_fields, n_format_fields, n_filter_fields);
	}

	/**<
	 *
	 * @param target
	 */
	inline void updateContainerSet(Index::BlockIndex::INDEX_BLOCK_TARGET target){
		// Determine target
		switch(target){
		case(Index::BlockIndex::INDEX_BLOCK_TARGET::INDEX_INFO)   :
			return(this->updateContainer(this->info_containers, this->index_entry.n_info_streams));
			break;
		case(Index::BlockIndex::INDEX_BLOCK_TARGET::INDEX_FORMAT) :
			return(this->updateContainer(this->format_containers, this->index_entry.n_format_streams));
			break;
		default: std::cerr << "unknown target type" << std::endl; exit(1);
		}
	}

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used word-size) for strides and data; if possible
	 *
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


	bool write(std::ofstream& stream, import_stats_type& stats, import_stats_type& stats_uncompressed);

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

	// Write everything
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.index_entry;

		if(entry.index_entry.controller.hasGTPermuted)
			stream << entry.ppa_manager;

		stream << entry.meta_hot_container;
		stream << entry.meta_cold_container;
		stream << entry.gt_rle_container;
		stream << entry.gt_simple_container;
		stream << entry.gt_support_data_container;
		stream << entry.meta_info_map_ids;
		stream << entry.meta_filter_map_ids;
		stream << entry.meta_format_map_ids;

		for(U32 i = 0; i < entry.index_entry.n_info_streams; ++i)
			stream << entry.info_containers[i];


		for(U32 i = 0; i < entry.index_entry.n_format_streams; ++i)
			stream << entry.format_containers[i];

		stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_BLOCK_EOF), sizeof(U64));

		return(stream);
	}

	// Read everything
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.index_entry;

		if(entry.index_entry.controller.hasGTPermuted) stream >> entry.ppa_manager;
		stream >> entry.meta_hot_container;
		stream >> entry.meta_cold_container;
		stream >> entry.gt_rle_container;
		stream >> entry.gt_simple_container;

		for(U32 i = 0; i < entry.index_entry.n_info_streams; ++i)
			stream >> entry.info_containers[i];


		for(U32 i = 0; i < entry.index_entry.n_format_streams; ++i)
			stream >> entry.format_containers[i];

		U64 eof_marker;
		stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
		assert(eof_marker == Constants::TACHYON_BLOCK_EOF);

		return(stream);
	}

public:
	index_entry_type  index_entry;
	permutation_type  ppa_manager;
	container_type    meta_hot_container;
	container_type    meta_info_map_ids;
	container_type    meta_format_map_ids;
	container_type    meta_filter_map_ids;
	container_type    meta_cold_container;
	container_type    gt_support_data_container; // data (0: rle, 1: simple), strides (n_objects)
	container_type    gt_rle_container;
	container_type    gt_simple_container;
	container_type*   info_containers;
	container_type*   format_containers;
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
