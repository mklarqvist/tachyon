#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include <unordered_map>

#include "third_party/xxhash/xxhash.h"

#include "components/variant_block_footer.h"
#include "components/variant_block_header.h"
#include "core/data_block_settings.h"
#include "data_container.h"
#include "core/meta_entry.h"
#include "core/variant_importer_container_stats.h"
#include "core/genotypes.h"

namespace tachyon{
namespace containers{

struct yon_blk_load_settings {
public:
	yon_blk_load_settings() : loaded_genotypes(false){}
	~yon_blk_load_settings(){}
	yon_blk_load_settings(const yon_blk_load_settings& other) :
		loaded_genotypes(other.loaded_genotypes),
		info_id_local_loaded(other.info_id_local_loaded),
		format_id_local_loaded(other.format_id_local_loaded),
		info_id_global_loaded(other.info_id_global_loaded),
		format_id_global_loaded(other.format_id_global_loaded),
		info_patterns_local(other.info_patterns_local),
		format_patterns_local(other.format_patterns_local),
		info_map_global(other.info_map_global),
		format_map_global(other.format_map_global)
	{

	}

	yon_blk_load_settings(yon_blk_load_settings&& other) noexcept :
		loaded_genotypes(other.loaded_genotypes),
		info_id_local_loaded(std::move(other.info_id_local_loaded)),
		format_id_local_loaded(std::move(other.format_id_local_loaded)),
		info_id_global_loaded(std::move(other.info_id_global_loaded)),
		format_id_global_loaded(std::move(other.format_id_global_loaded)),
		info_patterns_local(std::move(other.info_patterns_local)),
		format_patterns_local(std::move(other.format_patterns_local)),
		info_map_global(std::move(other.info_map_global)),
		format_map_global(std::move(other.format_map_global))
	{

	}

	yon_blk_load_settings& operator=(const yon_blk_load_settings& other){
		*this = yon_blk_load_settings(other);
		return(*this);
	}

	yon_blk_load_settings& operator=(yon_blk_load_settings&& other) noexcept{
		loaded_genotypes = other.loaded_genotypes;
		info_id_local_loaded = std::move(other.info_id_local_loaded);
		format_id_local_loaded = std::move(other.format_id_local_loaded);
		info_id_global_loaded = std::move(other.info_id_global_loaded);
		format_id_global_loaded = std::move(other.format_id_global_loaded);
		info_patterns_local = std::move(other.info_patterns_local);
		format_patterns_local = std::move(other.format_patterns_local);
		info_map_global = std::move(other.info_map_global);
		format_map_global = std::move(other.format_map_global);
		return(*this);
	}

	void clear(){
		this->loaded_genotypes = false;
		this->info_id_local_loaded.clear();
		this->format_id_local_loaded.clear();
		this->info_id_global_loaded.clear();
		this->format_id_global_loaded.clear();
		this->info_patterns_local.clear();
		this->format_patterns_local.clear();
		this->info_map_global.clear();
		this->format_map_global.clear();
	}

public:
	bool loaded_genotypes;
	std::vector<int> info_id_local_loaded;
	std::vector<int> format_id_local_loaded;
	std::vector<int> info_id_global_loaded;
	std::vector<int> format_id_global_loaded;
	std::vector< std::vector<int> > info_patterns_local;
	std::vector< std::vector<int> > format_patterns_local;
	std::unordered_map<int, int> info_map_global;
	std::unordered_map<int, int> format_map_global;
};

/**
 * Primary Tachyon block object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class VariantBlock {
public:
	typedef VariantBlock                    self_type;
	typedef DataContainer                   container_type;
	typedef VariantBlockHeader              block_header_type;
	typedef VariantBlockFooter              block_footer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef DataContainerHeader             offset_type;
	typedef tachyon::core::MetaEntry        meta_entry_type;
	typedef DataBlockSettings               block_settings_type;

public:
	VariantBlock();
	VariantBlock(const uint16_t n_info, const uint16_t n_format);
	~VariantBlock();
	VariantBlock(const self_type& other);
	VariantBlock(self_type&& other) noexcept;
	VariantBlock& operator=(const self_type& other);
	VariantBlock& operator=(self_type&& other) noexcept;

	void Allocate(const uint16_t n_info,
	              const uint16_t n_format,
	              const uint16_t n_filter);

	/**<
	 * Resize base container buffer streams
	 * @param s Number of bytes to allocate in buffers.
	 */
	void resize(const uint32_t s);

	/**<
	 * Recycle structure without releasing memory.
	 */
	void clear(void);

	inline const uint32_t& size(void) const{ return(this->header.n_variants); }

	/**<
	 * Reads all objects from disk. Primary function for reading
	 * entire blocks of data from disk. Data read in this way is
	 * not checked for integrity here.
	 * @param stream   Input stream
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool read(std::ifstream& stream);

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity here.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @param header   Reference global header.
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool read(std::ifstream& stream,
	          block_settings_type& settings,
	          const VariantHeader& header);

	/**<
	 * Read the header and footer of a block.
	 * @param stream
	 * @return
	 */
	bool ReadHeaderFooter(std::ifstream& stream);

	/**<
	 * Standard way of writing out a YON block.
	 * @param stream       Target output stream
	 * @param stats_basic  Tracking for basic containers
	 * @param stats_info   Tracking for INFO containers
	 * @param stats_format Tracking for FORMAT containers
	 * @return             Returns TRUE upon success or FALSE otherwise
	 */
	bool write(std::ostream& stream,
	           import_stats_type& stats_basic,
	           import_stats_type& stats_info,
	           import_stats_type& stats_format);

	/**<
	 * Add the data from a MetaEntry object to this block. Internally
	 * performs all the operations required to transfer each MetaEntry
	 * field into the correct destination with the correct encodings.
	 * This is the preferred way of storing a MetaEntry.
	 * @param meta_entry Input MetaEntry object to be stored.
	 * @return           Returns TRUE upon success or FALSE otherwise.
	 */
	bool operator+=(meta_entry_type& meta_entry);
	inline bool operator<<(meta_entry_type& meta_entry){ return(*this += meta_entry); }

	// Todo:
	bool AddInfo(const DataContainer& dc, container_type* dst_containers, int32_t(VariantBlock::*StreamFieldLookup)(int32_t)){
		if(dc.GetIdx() == -1){
			std::cerr << utility::timestamp("ERROR","IMPORT") << "DataContainer does not have not set the required field global idx." << std::endl;
			return false;
		}

		// 1) Search for global idx
		int32_t info_exists = (this->*StreamFieldLookup)(dc.GetIdx());
		if(info_exists >= 0){
			std::cerr << "field is already set. overwriting @ " << info_exists  << std::endl;
			dst_containers[info_exists] = dc;
		} else {
			std::cerr << "field does not exist. adding. " << info_exists  << std::endl;
		}

		return true;
	}

	/**<
	 * Compares a vector of global Info/Format/Filter identifiers to the identifier set in this
	 * block and returns the set intersection of keys.
	 * @param keys Vector of global Info/Format/Filter keys
	 * @return     Returns the set intersection of the provided keys and the local keys.
	 */
	std::vector<int> IntersectInfoKeys(const std::vector<int>& info_ids_global) const;
	std::vector<int> IntersectFormatKeys(const std::vector<int>& format_ids_global) const;
	std::vector<int> IntersectFilterKeys(const std::vector<int>& filter_ids_global) const;

	/**<
	 * Intersects a provided vector of global identifiers to a given pattern vector for
	 * a given Info/Format/Filter type.
	 * @param keys     Provided vector of global Info/Format/Filter keys.
	 * @param local_id Array offset to a local container.
	 * @return         Returns the set intersection of the provided keys and the target pattern keys.
	 */
	std::vector<int> IntersectInfoPatterns(const std::vector<int>& info_ids_global, const uint32_t local_id) const;
	std::vector<int> IntersectFormatPatterns(const std::vector<int>& format_ids_global, const uint32_t local_id) const;
	std::vector<int> IntersectFilterPatterns(const std::vector<int>& filter_ids_global, const uint32_t local_id) const;

	/**<
	 * Utility functions to retrieve all Info/Format/Filter keys from the footer as
	 * a vector of integers.
	 * @return Returns a vector of global identifiers.
	 */
	std::vector<uint32_t> GetInfoKeys(void) const;
	std::vector<uint32_t> GetFormatKeys(void) const;
	std::vector<uint32_t> GetFilterKeys(void) const;

	/**<
	 * Wrapper function to load a data container from packed YON blocks
	 * @param stream    Input file handler
	 * @param offset    Header object
	 * @param container Destination container object
	 * @return
	 */
	inline bool LoadContainer(std::ifstream& stream,
	                          const offset_type& offset,
	                          container_type& container)
	{
		container.header = offset;
		stream >> container;
		assert(container.header == offset);
		return(stream.good());
	}

	/**<
	 * Wrapper function to load a data container from packed YON blocks. Additionally
	 * performs a (potential) random seek to the start of the data sector before reading.
	 * @param stream    Input file handler
	 * @param offset    Header object
	 * @param container Destination container object
	 * @return
	 */
	inline bool LoadContainerSeek(std::ifstream& stream,
	                              const offset_type& offset,
	                              container_type& container)
	{
		stream.seekg(this->start_compressed_data_ + offset.data_header.offset);
		container.header = offset;
		stream >> container;
		assert(container.header == offset);
		return(stream.good());
	}

	/**<
	 * Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base container
	 * offsets and checks/builds.
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used primitive type) for strides and data; if possible
	 */
	void UpdateContainers(void);

	/**<
	 * Determine compressed block-size. Execute this function prior to writing a
	 * block.
	 * @return Returns the total number of bytes that will be written.
	 */
	uint64_t GetCompressedSize(void) const;

	inline void PackFooter(void){
		this->footer_support.reset();
		this->footer_support.buffer_data_uncompressed << this->footer;
		++this->footer_support;
	}

	inline uint32_t AddInfoPattern(const std::vector<int>& pattern){ return(this->footer.AddInfoPattern(pattern)); }
	inline uint32_t AddFormatPattern(const std::vector<int>& pattern){ return(this->footer.AddFormatPattern(pattern)); }
	inline uint32_t AddFilterPattern(const std::vector<int>& pattern){ return(this->footer.AddFilterPattern(pattern)); }
	inline uint32_t AddInfo(const uint32_t id){ return(this->footer.AddInfo(id)); }
	inline uint32_t AddFormat(const uint32_t id){ return(this->footer.AddFormat(id)); }
	inline uint32_t AddFilter(const uint32_t id){ return(this->footer.AddFilter(id)); }

	/**<
	 * Finalize this block for writing. Tells the header and footer objects to
	 * finish and then precomputes the virtual file offsets for each byte stream
	 * (container) and moves that offset to the footer. This operation is mandatory
	 * prior to writing this block!
	 */
	void Finalize(void){
		this->footer.Finalize();

		// Pre-calculate the virtual file offsets prior to writing the block
		// to disk. This is required during parallel processing as we do not
		// want to the writer slave to spend time compressing or doing any
		// other compute instead of simply writing at I/O saturated speeds.
		uint64_t b_offset = 0;
		if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
			this->UpdateHeader(this->footer.offsets[YON_BLK_PPA], this->base_containers[YON_BLK_PPA], 0);
			b_offset += this->base_containers[YON_BLK_PPA].GetObjectSize();
		}

		for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
			this->UpdateHeader(this->footer.offsets[i], this->base_containers[i], b_offset);
			b_offset += this->base_containers[i].GetObjectSize();
		}

		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
			this->UpdateHeader(this->footer.info_offsets[i], this->info_containers[i], b_offset);
			b_offset += this->info_containers[i].GetObjectSize();
		}

		for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
			this->UpdateHeader(this->footer.format_offsets[i], this->format_containers[i], b_offset);
			b_offset += this->format_containers[i].GetObjectSize();
		}
	}

	int32_t GetInfoPosition(const uint32_t global_id)   const;
	int32_t GetFormatPosition(const uint32_t global_id) const;
	int32_t GetFilterPosition(const uint32_t global_id) const;

	bool HasInfo(const uint32_t global_id)   const;
	bool HasFormat(const uint32_t global_id) const;
	bool HasFilter(const uint32_t global_id) const;

	container_type* GetInfoContainer(const uint32_t global_id)   const;
	container_type* GetFormatContainer(const uint32_t global_id) const;

	std::vector<bool> InfoPatternSetMembership(const int value)  const;
	std::vector<bool> FormatPatternSetMembership(const int value) const;
	std::vector<bool> FilterPatternSetMembership(const int value) const;

private:
	/**<
	 * Parse user-provided settings that provide information regarding
	 * how to slice and/or display data from this block prior to loading
	 * from disk/stream.
	 * @param settings Reference to user-provided settings.
	 * @param header   Reference to the global variant header object.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseSettings(DataBlockSettings& settings, const VariantHeader& header);

	/**<
	 * Parse what data will be displayed given the requested fields and
	 * what is available in the block. This is a private function as is
	 * called internally from ParseSettings().
	 * @param settings Reference to user-provided settings.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool ParseLoadedPatterns(DataBlockSettings& settings);

	/**<
	 *
	 * @param stats_basic
	 * @param stats_info
	 * @param stats_format
	 */
	void UpdateOutputStatistics(import_stats_type& stats_basic,
	                            import_stats_type& stats_info,
	                            import_stats_type& stats_format);

	/**<
	 * Move over pair of headers from a data container to a block footer
	 * @param offset    Destination header in footer
	 * @param container Target container hosting the header
	 */
	inline void UpdateHeader(offset_type& offset, const container_type& container){
		const uint32_t global_key = offset.data_header.global_key; // carry over global key
		offset = container.header;
		assert(offset == container.header); // Assert copy is correct
		offset.data_header.global_key = global_key;
	}

	/**<
	 * Move over pair of headers from a data container to a block footer
	 * @param offset         Destination header in footer
	 * @param container      Target container hosting the header
	 * @param virtual_offset Block virtual offset
	 */
	inline void UpdateHeader(offset_type& offset,
	                         const container_type& container,
	                         const uint32_t& virtual_offset)
	{
		const uint32_t global_key = offset.data_header.global_key; // carry over global key
		offset = container.header;
		assert(offset == container.header); // Assert copy is correct
		offset.data_header.global_key = global_key;
		offset.data_header.offset     = virtual_offset;
	}

	/**<
	 *
	 * @param stream
	 * @param offset
	 * @param container
	 * @param virtual_offset
	 */
	inline void WriteContainer(std::ostream& stream,
	                           offset_type& offset,
	                           const container_type& container,
	                           const uint32_t virtual_offset)
	{
		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
			return(this->WriteContainerEncrypted(stream, offset, container, virtual_offset));

		assert(offset.data_header.offset == virtual_offset);
		this->UpdateHeader(offset, container, virtual_offset);
		assert(container.buffer_data.size() == offset.data_header.cLength);
		stream << container;
	}

	/**<
	 *
	 * @param stream
	 * @param offset
	 * @param container
	 * @param virtual_offset
	 */
	inline void WriteContainerEncrypted(std::ostream& stream,
	                                    offset_type& offset,
	                                    const container_type& container,
	                                    const uint32_t virtual_offset)
	{
		this->UpdateHeader(offset, container, virtual_offset);
		assert(container.buffer_data.size() == offset.data_header.eLength);
		// Encrypted data is concatenated: write only data buffer
		stream.write(container.buffer_data.data(), container.buffer_data.size());
	}

public:
	uint16_t n_info_c_allocated;
	uint16_t n_format_c_allocated;
	block_header_type header;
	block_footer_type footer;
	container_type*   base_containers;
	container_type*   info_containers;
	container_type*   format_containers;
	yon_gt_ppa*       gt_ppa;
	yon_blk_load_settings* load_settings;

	// Utility
	uint64_t end_block_;
	uint64_t start_compressed_data_;
	uint64_t end_compressed_data_;
	container_type footer_support; // used internally only
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
