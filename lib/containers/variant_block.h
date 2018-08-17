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

/**
 * Primary Tachyon block object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class VariantBlock{
public:
	typedef VariantBlock                    self_type;
	typedef DataContainer                   container_type;
	typedef VariantBlockHeader              block_header_type;
	typedef VariantBlockFooter              block_footer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef DataContainerHeader             offset_type;
	typedef tachyon::core::MetaEntry        meta_entry_type;

public:
	VariantBlock();
	VariantBlock(const uint16_t n_info, const uint16_t n_format);
	~VariantBlock();

	VariantBlock(const self_type& other) :
		n_info_c_allocated(other.n_info_c_allocated),
		n_format_c_allocated(other.n_info_c_allocated),
		header(other.header),
		footer(other.footer),
		base_containers(new container_type[YON_BLK_N_STATIC]),
		info_containers(new container_type[other.n_info_c_allocated]),
		format_containers(new container_type[other.n_format_c_allocated]),
		gt_ppa(nullptr),
		end_block_(other.end_block_),
		start_compressed_data_(other.start_compressed_data_),
		end_compressed_data_(other.end_compressed_data_),
		footer_support(other.footer_support)
	{
		if(other.gt_ppa != nullptr){
			// Copy ppa data to new object.
			this->gt_ppa = new yon_gt_ppa(*other.gt_ppa);
		}

		for(U32 i = 0; i < YON_BLK_N_STATIC; ++i) this->base_containers[i] = other.base_containers[i];
		for(U32 i = 0; i < this->n_info_c_allocated; ++i)   this->info_containers[i]   = other.info_containers[i];
		for(U32 i = 0; i < this->n_format_c_allocated; ++i) this->format_containers[i] = other.format_containers[i];
	}

	VariantBlock(self_type&& other) noexcept :
		n_info_c_allocated(other.n_info_c_allocated),
		n_format_c_allocated(other.n_info_c_allocated),
		header(std::move(other.header)),
		footer(std::move(other.footer)),
		base_containers(nullptr),
		info_containers(nullptr),
		format_containers(nullptr),
		gt_ppa(nullptr),
		end_block_(other.end_block_),
		start_compressed_data_(other.start_compressed_data_),
		end_compressed_data_(other.end_compressed_data_),
		footer_support(std::move(other.footer_support))
	{
		std::swap(this->base_containers, other.base_containers);
		std::swap(this->info_containers, other.info_containers);
		std::swap(this->format_containers, other.format_containers);
		std::swap(this->gt_ppa, other.gt_ppa);
	}

	VariantBlock& operator=(const self_type& other){
		delete [] this->base_containers;
		delete [] this->info_containers;
		delete [] this->format_containers;
		delete this->gt_ppa;
		*this = VariantBlock(other);
		return(*this);
	}

	VariantBlock& operator=(self_type&& other) noexcept{
		if(this == &other){
			// precautions against self-moves
			return *this;
		}

		this->n_info_c_allocated = other.n_info_c_allocated;
		this->n_format_c_allocated = other.n_format_c_allocated;
		this->header = std::move(other.header);
		this->footer = std::move(other.footer);
		delete [] this->base_containers; this->base_containers = nullptr;
		std::swap(this->base_containers, other.base_containers);
		delete [] this->info_containers; this->info_containers = nullptr;
		std::swap(this->info_containers, other.info_containers);
		delete [] this->format_containers; this->format_containers = nullptr;
		std::swap(this->format_containers, other.format_containers);
		delete this->gt_ppa; this->gt_ppa = nullptr;
		std::swap(this->gt_ppa, other.gt_ppa);
		this->end_block_ = other.end_block_;
		this->start_compressed_data_ = other.start_compressed_data_;
		this->end_compressed_data_ = other.end_compressed_data_;
		this->footer_support = std::move(other.footer_support);
		return(*this);
	}

	void Allocate(const uint16_t n_info,
	              const uint16_t n_format,
	              const uint16_t n_filter)
	{
		// Allocate space for INFO containers.
		delete [] this->info_containers;
		this->info_containers = new container_type[n_info];
		this->n_info_c_allocated = n_info;

		// Allocate space for FORMAT containers.
		delete [] this->format_containers;
		this->format_containers = new container_type[n_format];
		this->n_format_c_allocated = n_format;

		// Alocate space for headers.
		this->footer.AllocateHeaders(n_info, n_format, n_filter);
	}

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

	/**< @brief Reads all digital objects from disk
	 * Primary function for reading data from disk. Data
	 * read in this way is not checked for integrity here.
	 * @param stream   Input stream
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool read(std::ifstream& stream);

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
	bool write(std::ostream& stream, import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format);

	// Add a meta entry
	/**<
	 * Overloaded operator for adding a variant entry
	 * @param meta_entry
	 * @return
	 */
	bool operator+=(meta_entry_type& meta_entry);
	inline bool operator<<(meta_entry_type& meta_entry){ return(*this += meta_entry); }

	/**<
	 * Compares a vector of global INFO identifiers to the identifier set in this
	 * block and returns the set intersection of keys
	 * @param info_ids Vector of global INFO keys
	 * @return         Returns the set intersection of provided keys and local keys
	 */
	std::vector<int> IntersectInfoKeys(const std::vector<int>& info_ids_global) const{
		std::vector<int> info_ids_found;
		if(info_ids_global.size() == 0) return(info_ids_found);

		for(U32 i = 0; i < info_ids_global.size(); ++i){
			for(U32 j = 0; j < this->footer.n_info_streams; ++j){
				if(this->footer.info_offsets[j].data_header.global_key == info_ids_global[i])
					info_ids_found.push_back(this->footer.info_offsets[j].data_header.global_key);
			}
		}

		return(info_ids_found);
	}

	std::vector<int> IntersectInfoPatterns(const std::vector<int>& info_ids_global, const uint32_t local_id) const{
		std::vector<int> info_ids_found;
		if(info_ids_global.size() == 0) return(info_ids_found);
		assert(local_id < this->footer.n_info_patterns);

		for(U32 i = 0; i < info_ids_global.size(); ++i){
			for(U32 k = 0; k < this->footer.info_patterns[local_id].pattern.size(); ++k){
				if(this->footer.info_patterns[local_id].pattern[k] == info_ids_global[i]){
					info_ids_found.push_back(this->footer.info_patterns[local_id].pattern[k]);
				}
			}
		}

		return(info_ids_found);
	}

	/**<
	 * Compares a vector of global FORMAT identifiers to the identifier set in this
	 * block and returns the set intersection of keys
	 * @param info_ids Vector of global FORMAT keys
	 * @return         Returns the set intersection of provided keys and local keys
	 */
	std::vector<int> IntersectFormatKeys(const std::vector<int>& format_ids_global) const{
		std::vector<int> format_ids_found;
		if(format_ids_global.size() == 0) return(format_ids_found);

		for(U32 i = 0; i < format_ids_global.size(); ++i){
			for(U32 j = 0; j < this->footer.n_format_streams; ++j){
				if(this->footer.format_offsets[j].data_header.global_key == format_ids_global[i])
					format_ids_found.push_back(this->footer.format_offsets[j].data_header.global_key);
			}
		}

		return(format_ids_found);
	}

	std::vector<int> IntersectFormatPatterns(const std::vector<int>& format_ids_global, const uint32_t local_id) const{
		std::vector<int> format_ids_found;
		if(format_ids_global.size() == 0) return(format_ids_found);
		assert(local_id < this->footer.n_format_patterns);

		for(U32 i = 0; i < format_ids_global.size(); ++i){
			for(U32 k = 0; k < this->footer.format_patterns[local_id].pattern.size(); ++k){
				if(this->footer.format_patterns[local_id].pattern[k] == format_ids_global[i])
					format_ids_found.push_back(this->footer.format_patterns[local_id].pattern[k]);
			}
		}

		return(format_ids_found);
	}

	std::vector<U32> GetFormatKeys(void) const{
		std::vector<U32> ret;
		for(U32 i = 0; i < this->footer.n_format_streams; ++i)
			ret.push_back(this->footer.format_offsets[i].data_header.global_key);

		return(ret);
	}

	std::vector<U32> GetInfoKeys(void) const{
		std::vector<U32> ret;
		for(U32 i = 0; i < this->footer.n_info_streams; ++i)
			ret.push_back(this->footer.info_offsets[i].data_header.global_key);

		return(ret);
	}

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

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used primitive type) for strides and data; if possible
	 */
	void UpdateContainers(void);

	/**<
	 * Determine compressed block-size. Execute this function prior to writing a
	 * block
	 * @return Returns the sum total disk size
	 */
	U64 DetermineCompressedSize(void) const;

	inline void PackFooter(void){
		this->footer_support.reset();
		this->footer_support.buffer_data_uncompressed << this->footer;
		++this->footer_support;
	}

	inline U32 AddInfoPattern(const std::vector<int>& pattern){ return(this->footer.AddInfoPattern(pattern)); }
	inline U32 AddFormatPattern(const std::vector<int>& pattern){ return(this->footer.AddFormatPattern(pattern)); }
	inline U32 AddFilterPattern(const std::vector<int>& pattern){ return(this->footer.AddFilterPattern(pattern)); }
	inline U32 AddInfo(const U32 id){ return(this->footer.AddInfo(id)); }
	inline U32 AddFormat(const U32 id){ return(this->footer.AddFormat(id)); }
	inline U32 AddFilter(const U32 id){ return(this->footer.AddFilter(id)); }
	inline void Finalize(void){ this->footer.Finalize(); }

	S32 GetInfoPosition(const U32 global_id) const{
		if(this->footer.info_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
		if(it == this->footer.info_map->end()) return -1;
		return(it->second);
	}

	S32 GetFormatPosition(const U32 global_id) const{
		if(this->footer.format_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
		if(it == this->footer.format_map->end()) return -1;
		return(it->second);
	}

	S32 GetFilterPosition(const U32 global_id) const{
		if(this->footer.filter_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
		if(it == this->footer.filter_map->end()) return -1;
		return(it->second);
	}

	bool HasInfo(const U32 global_id) const{
		if(this->footer.info_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
		if(it == this->footer.info_map->end()) return false;
		return(true);
	}

	bool HasFormat(const U32 global_id) const{
		if(this->footer.format_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
		if(it == this->footer.format_map->end()) return false;
		return(true);
	}

	bool HasFilter(const U32 global_id) const{
		if(this->footer.filter_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
		if(it == this->footer.filter_map->end()) return false;
		return(true);
	}

	container_type* GetInfoContainer(const U32 global_id) const{
		if(this->HasInfo(global_id))
			return(&this->info_containers[this->footer.info_map->at(global_id)]);
		else
			return nullptr;
	}

	container_type* GetFormatContainer(const U32 global_id) const{
		if(this->HasFormat(global_id))
			return(&this->format_containers[this->footer.format_map->at(global_id)]);
		else
			return nullptr;
	}

	std::vector<bool> InfoPatternSetMembership(const int value) const{
		std::vector<bool> matches(this->footer.n_info_patterns, false);
		for(U32 i = 0; i < this->footer.n_info_patterns; ++i){
			for(U32 j = 0; j < this->footer.info_patterns[i].pattern.size(); ++j){
				if(this->footer.info_patterns[i].pattern[j] == value){
					matches[i] = true;
					break;
				}
			}
		}
		return(matches);
	}

	std::vector<bool> FormatPatternSetMembership(const int value) const{
		std::vector<bool> matches(this->footer.n_format_patterns, false);
		for(U32 i = 0; i < this->footer.n_format_patterns; ++i){
			for(U32 j = 0; j < this->footer.format_patterns[i].pattern.size(); ++j){
				if(this->footer.format_patterns[i].pattern[j] == value){
					matches[i] = true;
					break;
				}
			}
		}
		return(matches);
	}

	std::vector<bool> FilterPatternSetMembership(const int value) const{
		std::vector<bool> matches(this->footer.n_filter_patterns, false);
		for(U32 i = 0; i < this->footer.n_filter_patterns; ++i){
			for(U32 j = 0; j < this->footer.filter_patterns[i].pattern.size(); ++j){
				if(this->footer.filter_patterns[i].pattern[j] == value){
					matches[i] = true;
					break;
				}
			}
		}
		return(matches);
	}

private:
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
		const U32 global_key = offset.data_header.global_key; // carry over global key
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
	                         const U32& virtual_offset)
	{
		const U32 global_key = offset.data_header.global_key; // carry over global key
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
	                           const U32 virtual_offset)
	{
		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
			return(this->WriteContainerEncrypted(stream, offset, container, virtual_offset));

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
	                                    const U32 virtual_offset)
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
	yon_gt_ppa* gt_ppa;

	// Utility
	uint64_t end_block_;
	uint64_t start_compressed_data_;
	uint64_t end_compressed_data_;
	container_type footer_support; // used internally only
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
