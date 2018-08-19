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

	bool readSlice(std::ifstream& stream){
		block_settings_type& settings;

		// Allocate memory for the Format and Info containers.
		// Info containers.
		delete [] this->info_containers;
		this->info_containers = new VariantBlock::container_type[this->footer.n_info_streams];
		this->n_info_c_allocated = this->footer.n_info_streams;
		// Format containers.
		delete [] this->format_containers;
		this->format_containers = new VariantBlock::container_type[this->footer.n_format_streams];
		this->n_format_c_allocated = this->footer.n_format_streams;

		// Interpret the user-specified block-settings if any. This step converts
		// global index offset values into local offsets and computes new pattern
		// vectors if required. The ordering of the values are according to the
		// input sequence not according to the actual stored order.
		this->ParseSettings(settings);

		// Load the FORMAT:GT (GBPBWT) permutation array.
		if(settings.load_static & YON_BLK_BV_PPA){
			// If there is FORMAT:GT field data available AND that data has
			// been permuted then create a new yon_gt_ppa object to store
			// this data.
			if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
				stream.seekg(this->start_compressed_data_ + this->footer.offsets[YON_BLK_PPA].data_header.offset);
				this->LoadContainerSeek(stream,
											   this->footer.offsets[YON_BLK_PPA],
											   this->base_containers[YON_BLK_PPA]);

				this->gt_ppa = new yon_gt_ppa;
				this->gt_ppa->n_s = this->header_->GetNumberSamples();
			}
		}

		// Load base meta containers.
		for(uint32_t i = YON_BLK_CONTIG; i < YON_BLK_GT_INT8; ++i){
			if(settings.load_static & (1 << i)){
				this->LoadContainerSeek(stream,
											   this->footer.offsets[i],
											   this->base_containers[i]);
			}
		}

		// Load genotype containers. At the moment, genotype containers cannot be loaded
		// individually by using this wrapper routine. If you wish to load these separately
		// you will have to do so manually.
		if((settings.load_static & YON_BLK_BV_GT) || (settings.load_static & YON_BLK_BV_FORMAT)){
			this->loaded_genotypes = true;
			this->LoadContainerSeek(stream, this->footer.offsets[YON_BLK_GT_INT8], this->base_containers[YON_BLK_GT_INT8]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT16],    this->base_containers[YON_BLK_GT_INT16]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT32],    this->base_containers[YON_BLK_GT_INT32]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT64],    this->base_containers[YON_BLK_GT_INT64]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT8],   this->base_containers[YON_BLK_GT_S_INT8]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT16],  this->base_containers[YON_BLK_GT_S_INT16]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT32],  this->base_containers[YON_BLK_GT_S_INT32]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT64],  this->base_containers[YON_BLK_GT_S_INT64]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT8],   this->base_containers[YON_BLK_GT_N_INT8]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT16],  this->base_containers[YON_BLK_GT_N_INT16]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT32],  this->base_containers[YON_BLK_GT_N_INT32]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT64],  this->base_containers[YON_BLK_GT_N_INT64]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_SUPPORT],  this->base_containers[YON_BLK_GT_SUPPORT]);
			this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_PLOIDY],   this->base_containers[YON_BLK_GT_PLOIDY]);
		}

		// Load Info containers. Technically there is no difference between the two
		// conditions below in terms of outcome. However, the first case guarantees
		// that data is loaded linearly from disk as this can be guaranteed when loading
		// all available data. There is no such guarntees for the second case.
		if(this->footer.n_info_streams && (settings.load_static & YON_BLK_BV_INFO) && settings.annotate_extra == false){
			stream.seekg(this->start_compressed_data_ + this->footer.info_offsets[0].data_header.offset);

			for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
				this->LoadContainer(stream,
										   this->footer.info_offsets[i],
										   this->info_containers[i]);
			}
		}
		// If we have a user-supplied list of identifiers parsed above.
		else {
			for(uint32_t i = 0; i < this->info_id_local_loaded.size(); ++i){
				this->LoadContainerSeek(stream,
				                        this->footer.info_offsets[this->info_id_local_loaded[i]],
				                        this->info_containers[this->info_id_local_loaded[i]]);
			}

		}

		// Load Format containers. Technically there is no difference between the two
		// conditions below in terms of outcome. However, the first case guarantees
		// that data is loaded linearly from disk as this can be guaranteed when loading
		// all available data. There is no such guarntees for the second case.
		if(this->footer.n_format_streams && (settings.load_static & YON_BLK_BV_FORMAT)){
			stream.seekg(this->start_compressed_data_ + this->footer.format_offsets[0].data_header.offset);
			for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
				this->LoadContainerSeek(stream, this->footer.format_offsets[i], this->format_containers[i]);
			}
			// At this point the stream should be located at the end-of-block
			// marker as the Format information is stored last.
			assert(this->end_compressed_data_ == (uint64_t)stream.tellg());
		}
		// If we have a user-supplied list of identifiers parsed above.
		else {
			for(uint32_t i = 0; i < this->format_id_local_loaded.size(); ++i){
				this->LoadContainerSeek(stream, this->footer.format_offsets[this->format_id_local_loaded[i]], this->format_containers[this->format_id_local_loaded[i]]);
			}
		}

		// Seek to end-of-block position.
		stream.seekg(this->end_block_);
		return(true);

	}

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

	/**<
	 * Compares a vector of global INFO identifiers to the identifier set in this
	 * block and returns the set intersection of keys
	 * @param info_ids Vector of global INFO keys
	 * @return         Returns the set intersection of provided keys and local keys
	 */
	std::vector<int> IntersectInfoKeys(const std::vector<int>& info_ids_global) const{
		std::vector<int> info_ids_found;
		if(info_ids_global.size() == 0) return(info_ids_found);

		for(uint32_t i = 0; i < info_ids_global.size(); ++i){
			for(uint32_t j = 0; j < this->footer.n_info_streams; ++j){
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

		for(uint32_t i = 0; i < info_ids_global.size(); ++i){
			for(uint32_t k = 0; k < this->footer.info_patterns[local_id].pattern.size(); ++k){
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

		for(uint32_t i = 0; i < format_ids_global.size(); ++i){
			for(uint32_t j = 0; j < this->footer.n_format_streams; ++j){
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

		for(uint32_t i = 0; i < format_ids_global.size(); ++i){
			for(uint32_t k = 0; k < this->footer.format_patterns[local_id].pattern.size(); ++k){
				if(this->footer.format_patterns[local_id].pattern[k] == format_ids_global[i])
					format_ids_found.push_back(this->footer.format_patterns[local_id].pattern[k]);
			}
		}

		return(format_ids_found);
	}

	std::vector<uint32_t> GetFormatKeys(void) const{
		std::vector<uint32_t> ret;
		for(uint32_t i = 0; i < this->footer.n_format_streams; ++i)
			ret.push_back(this->footer.format_offsets[i].data_header.global_key);

		return(ret);
	}

	std::vector<uint32_t> GetInfoKeys(void) const{
		std::vector<uint32_t> ret;
		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)
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
	uint64_t DetermineCompressedSize(void) const;

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
	inline void Finalize(void){ this->footer.Finalize(); }

	int32_t GetInfoPosition(const uint32_t global_id) const{
		if(this->footer.info_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
		if(it == this->footer.info_map->end()) return -1;
		return(it->second);
	}

	int32_t GetFormatPosition(const uint32_t global_id) const{
		if(this->footer.format_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
		if(it == this->footer.format_map->end()) return -1;
		return(it->second);
	}

	int32_t GetFilterPosition(const uint32_t global_id) const{
		if(this->footer.filter_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
		if(it == this->footer.filter_map->end()) return -1;
		return(it->second);
	}

	bool HasInfo(const uint32_t global_id) const{
		if(this->footer.info_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
		if(it == this->footer.info_map->end()) return false;
		return(true);
	}

	bool HasFormat(const uint32_t global_id) const{
		if(this->footer.format_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
		if(it == this->footer.format_map->end()) return false;
		return(true);
	}

	bool HasFilter(const uint32_t global_id) const{
		if(this->footer.filter_map == nullptr) return false;
		VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
		if(it == this->footer.filter_map->end()) return false;
		return(true);
	}

	container_type* GetInfoContainer(const uint32_t global_id) const{
		if(this->HasInfo(global_id))
			return(&this->info_containers[this->footer.info_map->at(global_id)]);
		else
			return nullptr;
	}

	container_type* GetFormatContainer(const uint32_t global_id) const{
		if(this->HasFormat(global_id))
			return(&this->format_containers[this->footer.format_map->at(global_id)]);
		else
			return nullptr;
	}

	std::vector<bool> InfoPatternSetMembership(const int value) const{
		std::vector<bool> matches(this->footer.n_info_patterns, false);
		for(uint32_t i = 0; i < this->footer.n_info_patterns; ++i){
			for(uint32_t j = 0; j < this->footer.info_patterns[i].pattern.size(); ++j){
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
		for(uint32_t i = 0; i < this->footer.n_format_patterns; ++i){
			for(uint32_t j = 0; j < this->footer.format_patterns[i].pattern.size(); ++j){
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
		for(uint32_t i = 0; i < this->footer.n_filter_patterns; ++i){
			for(uint32_t j = 0; j < this->footer.filter_patterns[i].pattern.size(); ++j){
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
