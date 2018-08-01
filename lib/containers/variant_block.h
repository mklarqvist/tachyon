#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include <unordered_map>

#include "third_party/xxhash/xxhash.h"

#include "algorithm/permutation/permutation_manager.h"
#include "components/variant_block_footer.h"
#include "components/variant_block_header.h"
#include "core/data_block_settings.h"
#include "data_container.h"
#include "core/meta_entry.h"
#include "core/variant_importer_container_stats.h"

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
	typedef VariantBlockHeader              block_header_type;
	typedef VariantBlockFooter              block_footer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef DataContainerHeader             offset_type;
	typedef tachyon::core::MetaEntry        meta_entry_type;

public:
	VariantBlock();
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
	std::vector<U32> IntersectInfoKeys(const std::vector<U32>& info_ids) const{
		std::vector<U32> info_ids_found;
		if(info_ids.size() == 0) return(info_ids_found);

		for(U32 i = 0; i < info_ids.size(); ++i){
			for(U32 j = 0; j < this->footer.n_info_streams; ++j){
				if(this->footer.info_offsets[j].data_header.global_key == info_ids[i])
					info_ids_found.push_back(this->footer.info_offsets[j].data_header.global_key);
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
	std::vector<U32> IntersectFormatKeys(const std::vector<U32>& format_ids) const{
		std::vector<U32> format_ids_found;
		if(format_ids.size() == 0) return(format_ids_found);

		for(U32 i = 0; i < format_ids.size(); ++i){
			for(U32 j = 0; j < this->footer.n_format_streams; ++j){
				if(this->footer.format_offsets[j].data_header.global_key == format_ids[i])
					format_ids_found.push_back(this->footer.format_offsets[j].data_header.global_key);
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
	inline bool LoadContainer(std::ifstream& stream, const offset_type& offset, container_type& container){
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
	inline bool LoadContainerSeek(std::ifstream& stream, const offset_type& offset, container_type& container){
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
	void UpdateOutputStatistics(import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format);

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
	inline void UpdateHeader(offset_type& offset, const container_type& container, const U32& virtual_offset){
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
	inline void WriteContainer(std::ostream& stream, offset_type& offset, const container_type& container, const U32 virtual_offset){
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
	inline void WriteContainerEncrypted(std::ostream& stream, offset_type& offset, const container_type& container, const U32 virtual_offset){
		this->UpdateHeader(offset, container, virtual_offset);
		assert(container.buffer_data.size() == offset.data_header.eLength);
		// Encrypted data is concatenated: write only data buffer
		stream.write(container.buffer_data.data(), container.buffer_data.size());
	}

public:
	block_header_type header;
	block_footer_type footer;
	permutation_type  ppa_manager; // Todo: exchange for yon_gt_ppa
	container_type*   base_containers;
	container_type*   info_containers;
	container_type*   format_containers;

	// Utility
	U64 end_block_;
	U64 start_compressed_data_;
	U64 end_compressed_data_;
	container_type footer_support; // used internally only
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
