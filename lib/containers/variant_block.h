#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include "algorithm/permutation/permutation_manager.h"
#include "components/variant_block_footer.h"
#include "components/variant_block_header.h"
#include "core/data_block_settings.h"
#include "data_container.h"
#include "core/meta_entry.h"
#include "core/variant_importer_container_stats.h"
#include "io/vcf/VCFHeader.h"

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
	typedef HashContainer                   hash_container_type;
	typedef HashVectorContainer             hash_vector_container_type;
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

	//
	//inline const size_t getINFOLoaded(void) const{ return(this->info_loaded.size()); }
	//inline const size_t getFORMATLoaded(void) const{ return(this->format_loaded.size()); }

	inline U32 addFieldINFO(const U32 fieldID){ return(this->info_fields.setGet(fieldID)); }
	inline U32 addFieldFORMAT(const U32 fieldID){ return(this->format_fields.setGet(fieldID)); }
	inline U32 addFieldFILTER(const U32 fieldID){ return(this->filter_fields.setGet(fieldID)); }

	inline const S32 getPatternsINFO(const U64& hash_pattern) const{
		U32 mapID = 0;
		if(this->info_patterns.getRaw(hash_pattern, mapID))
			return(mapID);
		else return(-1);
	}

	inline const S32 getPatternsFORMAT(const U64& hash_pattern) const{
		U32 mapID = 0;
		if(this->format_patterns.getRaw(hash_pattern, mapID))
			return(mapID);
		else return(-1);
	}

	inline const S32 getPatternsFILTER(const U64& hash_pattern) const{
		U32 mapID = 0;
		if(this->filter_patterns.getRaw(hash_pattern, mapID))
			return(mapID);
		else return(-1);
	}

	inline void addPatternINFO(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->info_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			exit(1);
		}
	}

	inline void addPatternFORMAT(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->format_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			exit(1);
		}
	}

	inline void addPatternFILTER(const std::vector<U32>& pattern, const U64& hash_pattern){
		if(!this->filter_patterns.set(pattern, hash_pattern)){
			std::cerr << "failed to insert filter: " << pattern.size() << " and " << hash_pattern << std::endl;
			exit(1);
		}
	}

	/**<
	 * Finalize this block before writing to disk. This wrapper function
	 * calls all necessary functions to construct a valid Tachyon block
	 * for sequence variant data
	 */
	inline void finalize(void){
		this->footer.n_info_streams   = this->info_fields.size();
		this->footer.n_filter_streams = this->filter_fields.size();
		this->footer.n_format_streams = this->format_fields.size();
		this->footer.allocateDiskOffsets(this->footer.n_info_streams, this->footer.n_format_streams, this->footer.n_filter_streams);
		this->updateContainers();
		this->footer.constructBitVector(containers::VariantBlockFooter::INDEX_INFO,   this->info_fields,   this->info_patterns);
		this->footer.constructBitVector(containers::VariantBlockFooter::INDEX_FILTER, this->filter_fields, this->filter_patterns);
		this->footer.constructBitVector(containers::VariantBlockFooter::INDEX_FORMAT, this->format_fields, this->format_patterns);
	}

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
	bool readHeaderFooter(std::ifstream& stream);

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
	std::vector<U32> intersectInfoKeys(const std::vector<U32>& info_ids) const{
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
	std::vector<U32> intersectFormatKeys(const std::vector<U32>& format_ids) const{
		std::vector<U32> format_ids_found;
		if(format_ids.size() == 0) return(format_ids_found);

		for(U32 i = 0; i < format_ids.size(); ++i){
			for(U32 j = 0; j < this->footer.n_info_streams; ++j){
				if(this->footer.info_offsets[j].data_header.global_key == format_ids[i])
					format_ids_found.push_back(this->footer.info_offsets[j].data_header.global_key);
			}
		}

		return(format_ids_found);
	}

	std::vector<U32> getFormatKeys(void) const{
		std::vector<U32> ret;
		for(U32 i = 0; i < this->footer.n_format_streams; ++i)
			ret.push_back(this->footer.format_offsets[i].data_header.global_key);

		return(ret);
	}

	std::vector<U32> getInfoKeys(void) const{
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
	inline bool __loadContainer(std::ifstream& stream, const offset_type& offset, container_type& container){
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
	inline bool __loadContainerSeek(std::ifstream& stream, const offset_type& offset, container_type& container){
		stream.seekg(this->start_compressed_data_ + offset.data_header.offset);
		container.header = offset;
		stream >> container;
		assert(container.header == offset);
		return(stream.good());
	}

private:
	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used primitive type) for strides and data; if possible
	 */
	void updateContainers(void);

	/**<
	 * Determine compressed block-size. Execute this function prior to writing a
	 * block
	 * @return Returns the sum total disk size
	 */
	const U64 __determineCompressedSize(void) const;

	/**<
	 *
	 * @param stats_basic
	 * @param stats_info
	 * @param stats_format
	 */
	void updateOutputStatistics(import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format);

	/**<
	 * Move over pair of headers from a data container to a block footer
	 * @param offset    Destination header in footer
	 * @param container Target container hosting the header
	 */
	inline void __updateHeader(offset_type& offset, const container_type& container){
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
	inline void __updateHeader(offset_type& offset, const container_type& container, const U32& virtual_offset){
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
	inline void __writeContainer(std::ostream& stream, offset_type& offset, const container_type& container, const U32 virtual_offset){
		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
			return(this->__writeContainerEncrypted(stream, offset, container, virtual_offset));

		this->__updateHeader(offset, container, virtual_offset);
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
	inline void __writeContainerEncrypted(std::ostream& stream, offset_type& offset, const container_type& container, const U32 virtual_offset){
		this->__updateHeader(offset, container, virtual_offset);
		assert(container.buffer_data.size() == offset.data_header.eLength);
		// Encrypted data is concatenated: write only data buffer
		stream.write(container.buffer_data.data(), container.buffer_data.size());
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
	container_type    gt_support_data_container;
	container_type    gt_rle8_container;
	container_type    gt_rle16_container;
	container_type    gt_rle32_container;
	container_type    gt_rle64_container;
	container_type    gt_simple8_container;
	container_type    gt_simple16_container;
	container_type    gt_simple32_container;
	container_type    gt_simple64_container;
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
	U64 end_block_;
	U64 start_compressed_data_;
	U64 end_compressed_data_;
	container_type footer_support; // used internally only
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
