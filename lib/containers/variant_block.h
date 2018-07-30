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
#include "io/vcf/vcf_header.h"

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
	typedef std::unordered_map<U32, U32>    map_type;
	typedef std::unordered_map<U64, U32>    map_pattern_type;

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

	inline U32 addFieldINFO(const U32 fieldID){ return(this->info_fields.setGet(fieldID)); }
	inline U32 addFieldFORMAT(const U32 fieldID){ return(this->format_fields.setGet(fieldID)); }
	inline U32 addFieldFILTER(const U32 fieldID){ return(this->filter_fields.setGet(fieldID)); }

	inline S32 getPatternsINFO(const U64& hash_pattern) const{
		U32 mapID = 0;
		if(this->info_patterns.getRaw(hash_pattern, mapID))
			return(mapID);
		else return(-1);
	}

	inline S32 getPatternsFORMAT(const U64& hash_pattern) const{
		U32 mapID = 0;
		if(this->format_patterns.getRaw(hash_pattern, mapID))
			return(mapID);
		else return(-1);
	}

	inline S32 getPatternsFILTER(const U64& hash_pattern) const{
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
	 * Todo: delete
	 * Finalize this block before writing to disk. This wrapper function
	 * calls all necessary functions to construct a valid Tachyon block
	 * for sequence variant data
	 */
	inline void finalize(void){
		this->footer.n_info_streams   = this->info_fields.size();
		this->footer.n_filter_streams = this->filter_fields.size();
		this->footer.n_format_streams = this->format_fields.size();
		this->footer.AllocateHeaders(this->footer.n_info_streams, this->footer.n_format_streams, this->footer.n_filter_streams);
		this->UpdateContainers();
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
			for(U32 j = 0; j < this->footer.n_format_streams; ++j){
				if(this->footer.format_offsets[j].data_header.global_key == format_ids[i])
					format_ids_found.push_back(this->footer.format_offsets[j].data_header.global_key);
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

	void Finalize(void){
		this->footer.ConstructInfoBitVector(this->info_map);
		this->footer.ConstructFormatBitVector(this->format_map);
		this->footer.ConstructFilterBitVector(this->filter_map);
	}

	bool BuildMaps(void){
		delete this->info_map;
		delete this->filter_map;
		delete this->format_map;

		this->info_map   = new map_type();
		this->filter_map = new map_type();
		this->format_map = new map_type();

		return true;
	}

	bool BuildPatternMaps(void){
		delete this->info_pattern_map;
		this->info_pattern_map = new map_pattern_type();
		if(this->footer.n_info_patterns_allocated == 0){
			this->footer.info_patterns = new yon_blk_bv_pair[100];
			this->footer.n_info_patterns_allocated = 100;
		}

		delete this->filter_pattern_map;
		this->filter_pattern_map = new map_pattern_type();
		if(this->footer.n_filter_patterns_allocated == 0){
			this->footer.filter_patterns = new yon_blk_bv_pair[100];
			this->footer.n_filter_patterns_allocated = 100;
		}

		delete this->format_pattern_map;
		this->format_pattern_map = new map_pattern_type();
		if(this->footer.n_format_patterns_allocated == 0){
			this->footer.format_patterns = new yon_blk_bv_pair[100];
			this->footer.n_format_patterns_allocated = 100;
		}

		return true;
	}

	U32 AddStreamWrapper(const U32 id, map_type* map, offset_type*& offsets, U16& stream_counter){
		map_type::const_iterator it = map->find(id);
		if(it == map->end()){
			(*map)[id] = stream_counter;
			offsets[stream_counter].data_header.global_key = id;
			++stream_counter;
		}

		return((*map)[id]);
	}

	U32 AddInfo(const U32 id){
		if(this->info_map == nullptr) this->BuildMaps();
		return(this->AddStreamWrapper(id, this->info_map, this->footer.info_offsets, this->footer.n_info_streams));
	}

	U32 AddFormat(const U32 id){
		if(this->format_map == nullptr) this->BuildMaps();
		return(this->AddStreamWrapper(id, this->format_map, this->footer.format_offsets, this->footer.n_format_streams));
	}

	U32 AddFilter(const U32 id){
		if(this->filter_map == nullptr) this->BuildMaps();
		return(this->AddStreamWrapper(id, this->filter_map, this->footer.filter_offsets, this->footer.n_filter_streams));
	}

	/**<
	* Static function that calculates the 64-bit hash value for the target
	* FORMAT/FILTER/INFO vector of id fields. The id fields must be of type
	* int (S32). Example of using this function:
	*
	* const U64 hash_value = VariantImporter::HashIdentifiers(id_vector);
	*
	* @param id_vector Input vector of FORMAT/FILTER/INFO identifiers.
	* @return          Returns a 64-bit hash value.
	*/
	static U64 HashIdentifiers(const std::vector<int>& id_vector){
		XXH64_state_t* const state = XXH64_createState();
		if (state==NULL) abort();

		XXH_errorcode const resetResult = XXH64_reset(state, BCF_HASH_SEED);
		if (resetResult == XXH_ERROR) abort();

		for(U32 i = 0; i < id_vector.size(); ++i){
			XXH_errorcode const addResult = XXH64_update(state, (const void*)&id_vector[i], sizeof(int));
			if (addResult == XXH_ERROR) abort();
		}

		U64 hash = XXH64_digest(state);
		XXH64_freeState(state);

		return hash;
	}

	U32 AddPatternWrapper(const std::vector<int>& pattern, map_pattern_type* pattern_map, yon_blk_bv_pair* bv_pairs, U16& stream_counter){
		U64 pattern_hash = VariantBlock::HashIdentifiers(pattern);
		const map_pattern_type::const_iterator it = pattern_map->find(pattern_hash); // search for pattern
		if(it == pattern_map->end()){
			(*pattern_map)[pattern_hash] = stream_counter;
			bv_pairs[stream_counter].pattern = pattern;
			++stream_counter;
		}

		return((*pattern_map)[pattern_hash]);
	}

	U32 AddInfoPattern(const std::vector<int>& pattern){
		if(this->info_pattern_map == nullptr) this->BuildPatternMaps();
		if(this->footer.n_info_patterns_allocated == 0){
			this->footer.info_patterns = new yon_blk_bv_pair[100];
			this->footer.n_info_patterns_allocated = 100;
		}
		return(this->AddPatternWrapper(pattern, this->info_pattern_map, this->footer.info_patterns, this->footer.n_info_patterns));
	}

	U32 AddFormatPattern(const std::vector<int>& pattern){
		if(this->format_pattern_map == nullptr) this->BuildPatternMaps();
		if(this->footer.n_format_patterns_allocated == 0){
			this->footer.format_patterns = new yon_blk_bv_pair[100];
			this->footer.n_format_patterns_allocated = 100;
		}
		return(this->AddPatternWrapper(pattern, this->format_pattern_map, this->footer.format_patterns, this->footer.n_format_patterns));
	}

	U32 AddFilterPattern(const std::vector<int>& pattern){
		if(this->filter_pattern_map == nullptr) this->BuildPatternMaps();
		if(this->footer.n_filter_patterns_allocated == 0){
			this->footer.filter_patterns = new yon_blk_bv_pair[100];
			this->footer.n_filter_patterns_allocated = 100;
		}
		return(this->AddPatternWrapper(pattern, this->filter_pattern_map, this->footer.filter_patterns, this->footer.n_filter_patterns));
	}

	/**<
	 * Determine compressed block-size. Execute this function prior to writing a
	 * block
	 * @return Returns the sum total disk size
	 */
	U64 DetermineCompressedSize(void) const;

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
			return(this->__writeContainerEncrypted(stream, offset, container, virtual_offset));

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
	inline void __writeContainerEncrypted(std::ostream& stream, offset_type& offset, const container_type& container, const U32 virtual_offset){
		this->UpdateHeader(offset, container, virtual_offset);
		assert(container.buffer_data.size() == offset.data_header.eLength);
		// Encrypted data is concatenated: write only data buffer
		stream.write(container.buffer_data.data(), container.buffer_data.size());
	}

public:
	block_header_type header;
	block_footer_type footer;
	permutation_type  ppa_manager;
	container_type*   base_containers;
	container_type*   info_containers;
	container_type*   format_containers;

	// Supportive hash tables to permit the map from global
	// IDX fields to local IDX fields.
	map_type* info_map;
	map_type* format_map;
	map_type* filter_map;
	map_pattern_type* info_pattern_map;
	map_pattern_type* format_pattern_map;
	map_pattern_type* filter_pattern_map;

	// Use during construction
	// Todo: Delete these
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
