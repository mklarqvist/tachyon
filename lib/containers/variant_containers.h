#ifndef CONTAINERS_VARIANT_CONTAINERS_H_
#define CONTAINERS_VARIANT_CONTAINERS_H_

#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <unordered_map>
#include <cmath>
#include <fstream>

#include "index/index.h"
#include "data_containers.h"
#include "primitive_container.h"
#include "core/variant_record.h"

namespace tachyon {

/** @brief Controller flags for an IndexBlockEntry
 * This structure is for internal use only and describes
 * various internal states as flags.
 */
struct yon_vb_hdr_cont {
	typedef yon_vb_hdr_cont self_type;

public:
	yon_vb_hdr_cont();
	~yon_vb_hdr_cont();

	void clear();

	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller);
	friend std::istream& operator>>(std::istream& stream, self_type& controller);

public:
	uint16_t
		has_gt:           1,  // This block has GT FORMAT data
		has_gt_permuted:  1,  // have the GT fields been permuted
		any_encrypted:    1,  // any data encrypted
		unused:          13; // reserved for future use
};

/** @brief Fixed-sized components of an IndexBlockEntry
 * For internal use only. This is a subcomponent of fixed
 * size describing:
 * 1) Contig, minimum and maximum genomic coordinates
 * 2) Offsets into the containers
 * 3) Number of containers and ID patterns
 * 4) Controller flags
 */
struct yon_vb_hdr{
private:
	typedef yon_vb_hdr           self_type;
	typedef yon_vb_hdr_cont controller_type;
	typedef yon_dc_hdr          header_type;

public:
	yon_vb_hdr();
	~yon_vb_hdr();

	inline const uint32_t& size(void) const{ return(this->n_variants); }
	inline const int32_t& GetContigId(void) const{ return(this->contig_id); }
	inline const int64_t& GetMinPosition(void) const{ return(this->min_position); }
	inline const int64_t& GetMaxPosition(void) const{ return(this->max_position); }
	inline uint64_t& GetBlockHash(void){ return(this->block_hash); }
	inline const uint64_t& GetBlockHash(void) const{ return(this->block_hash); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry);

	void reset(void);

public:
	// Allows jumping to the next block when streaming.
	// EOF marker is at this position - sizeof(EOF marker)
	uint32_t l_offset_footer;
	uint64_t block_hash;     // block identifier in the form of a random hash
	controller_type controller;
	int32_t  contig_id;       // contig identifier
	int64_t  min_position;    // minimum coordinate in this block
	int64_t  max_position;    // maximum coordinate in this block
	uint32_t n_variants;    // number of variants in this block
};

struct yon_vb_ftr {
public:
	typedef yon_vb_ftr  self_type;
	typedef yon_dc_hdr header_type;
	typedef std::unordered_map<uint32_t, uint32_t> map_type;
	typedef std::unordered_map<uint64_t, uint32_t> map_pattern_type;

public:
	yon_vb_ftr();
	~yon_vb_ftr();
	yon_vb_ftr(const self_type& other);
	yon_vb_ftr(self_type&& other) noexcept;
	yon_vb_ftr& operator=(const self_type& other);
	yon_vb_ftr& operator=(self_type&& other) noexcept;

	void reset(void);
	void resetTables(void);

	/**<
	 * Wrapper function for allocating memory for new offset objects
	 * for Info, Format, and Filter patterns and streams
	 * @param n_info_streams   Number of unique info streams.
	 * @param n_format_streams Number of unique format streams.
	 * @param n_filter_streams Number of unique filter streams.
	 */
	void AllocateHeaders(const uint32_t n_info_streams,
		                 const uint32_t n_format_streams,
		                 const uint32_t n_filter_streams);

	void AllocateInfoHeaders(const uint32_t n_info_streams);
	void AllocateFormatHeaders(const uint32_t n_format_streams);
	void AllocateFilterHeaders(const uint32_t n_filter_streams);

	/**<
	 * Wrapper function for adding a pattern to the block. The
	 * integers in the input vector is first hashed and check against
	 * the map of existing patterns. If the pattern does not exist
	 * then add it and return the local idx. Otherwise return the
	 * local idx for this pattern.
	 *
	 * Do not directly call this wrapper.
	 *
	 * Invocation of this function comes from the functions
	 * AddInfoPattern(), AddFormatPattern(), and
	 * AddFilterPattern().
	 *
	 * @param pattern        Input vector of global idx values.
	 * @param pattern_map    Target reference map of hashed global idx values.
	 * @param bv_pairs       Dst pointer of bit-vector entries.
	 * @param stream_counter Reference of current number of unique hash patterns.
	 * @return               Returns an array offset to the matching pattern (could be newly created).
	 */
	uint32_t AddPatternWrapper(const std::vector<int>& pattern,
	                           map_pattern_type* pattern_map,
	                           yon_blk_bv_pair* bv_pairs,
	                           uint16_t& stream_counter);
	uint32_t AddInfoPattern(const std::vector<int>& pattern);
	uint32_t AddFormatPattern(const std::vector<int>& pattern);
	uint32_t AddFilterPattern(const std::vector<int>& pattern);

	/**<
	 * Wrapper function to add a new byte stream to the block. Takes
	 * a input a global idx for a field and checks if that field is
	 * already set. If it is not then set it at the next available
	 * position and return that local idx. Otherwise, return the local
	 * idx of where this global idx has been set.
	 *
	 * Do not directly call this wrapper.
	 *
	 * Invocation of this function comes from the functions
	 * AddInfo(), AddFormat(), and
	 * AddFilter().
	 *
	 * @param global_id Input global idx.
	 * @param map       Target map from global to local idx.
	 * @param offsets   Pointer to dst offsets.
	 * @param n_streams Reference to number of currently set fields.
	 * @return          Returns the local idx for this global idx.
	 */
	uint32_t AddStreamWrapper(const uint32_t global_id,
	                          map_type* map,
	                          header_type*& offsets,
	                          uint16_t& n_streams);
	uint32_t AddInfo(const uint32_t global_id);
	uint32_t AddFormat(const uint32_t global_id);
	uint32_t AddFilter(const uint32_t global_id);

	/**<
	 * Perform operations required prior to writing disk to an output stream.
	 * At the moment, this involves constructing bit-vectors.
	 */
	void Finalize(void);

	/**<
	* Static function that calculates the 64-bit hash value for the target
	* FORMAT/FILTER/INFO vector of id fields. The id fields must be of type
	* int (S32). Example of using this function:
	*
	* const uint64_t hash_value = VariantImporter::HashIdentifiers(id_vector);
	*
	* @param id_vector Input vector of FORMAT/FILTER/INFO identifiers.
	* @return          Returns a 64-bit hash value.
	*/
	static uint64_t HashIdentifiers(const std::vector<int>& id_vector);

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& entry);

private:
	/**<
	 * Allocation functions for creating the stream and pattern maps
	 * from global idx to local idx.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool AllocateMaps(void);
	bool AllocatePatternMaps(void);

	/**<
	 * This wrapper adds patterns to the hash map when the data has
	 * already been loaded. This occurs when loading an object from
	 * disk/buffer.
	 *
	 * Do not directly call this wrapper.
	 *
	 * Invocation of this function comes from the functions
	 * UpdateInfoPatternMap(), UpdateFormatPatternMap(), and
	 * UpdateFilterPatternMap().
	 *
	 * @param pattern        Input vector of global idxs for the target field.
	 * @param pattern_map    Pointer to target map from global idx to local idx.
	 * @param stream_counter Reference target number of patterns to expect.
	 * @return               Returns the target local idx position where the pattern was set.
	 */
	uint32_t UpdatePatternMapWrapper(const std::vector<int>& pattern,
								     map_pattern_type* pattern_map,
								     const uint16_t& local_position);

	uint32_t UpdateInfoPatternMap(const std::vector<int>& pattern, const uint16_t local_position);
	uint32_t UpdateFormatPatternMap(const std::vector<int>& pattern, const uint16_t local_position);
	uint32_t UpdateFilterPatternMap(const std::vector<int>& pattern, const uint16_t local_position);

	/**<
	 * Update a given stream offset in the situation where the streams
	 * have already been loaded (during IO operations) and the number
	 * of fields are already known. The purpose of this function is to
	 * iterate over those fields and update the map from global to local
	 * idx.
	 *
	 * Do not directly call this wrapper.
	 *
	 * Invocation of this function comes from the functions
	 * UpdateInfoMap(), UpdateFormatMap(), and
	 * UpdateFilterMap().
	 *
	 * @param offset         Source header that has been preloaded and contains valid global idx information.
	 * @param map            Pointer to target map from global to local idx.
	 * @param local_position Local idx (array offset).
	 * @return               Returns the local idx offset.
	 */
	uint32_t UpdateOffsetMapWrapper(const header_type& offset,
									map_type* map,
									const uint16_t& local_position);
	uint32_t UpdateInfoMap(const header_type& offset, const uint16_t local_position);
	uint32_t UpdateFormatMap(const header_type& offset, const uint16_t local_position);
	uint32_t UpdateFilterMap(const header_type& offset, const uint16_t local_position);

	/**<
	 * Wrappers for constructing new Info/Format/Filter bit-vectors
	 * in the footer. This function requires you to pass the map
	 * from global idx values to local idx values. These wrappers
	 * are called from the Finalize() function when finishing a block
	 * for writing.
	 * @param pattern_map Pointer to target map from global idx to local idx.
	 * @return            Returns TRUE upon success or FALSE otherwise.
	 */
	bool ConstructInfoBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map);
	bool ConstructFormatBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map);
	bool ConstructFilterBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map);

public:
	// Utility members. l_*_bitvector stores the uint8_t-length
	// of each of the bit-vectors described below.
	// Not written or read from disk. Used internally only
	uint8_t l_info_bitvector;
	uint8_t l_format_bitvector;
	uint8_t l_filter_bitvector;

	// Critical values used to track the number of data streams
	// that is available for each possible type. The n_*_streams
	// fields corresponds to the number of different uint8_t streams
	// that are set in this block. The n_*_patterns corresponds
	// to the number of unique vectors of field identifiers that
	// occurred in the block.
	uint16_t n_info_streams; // streams
	uint16_t n_format_streams;
	uint16_t n_filter_streams;
	uint16_t n_info_patterns; // patterns
	uint16_t n_format_patterns;
	uint16_t n_filter_patterns;

	// Header structures corresponds critical information regarding
	// the global IDX and virtual uint8_t offset to the start of each
	// compressed and possibly encrypted uint8_t stream. In addition,
	// this structure details the primitive type of the data in the
	// stream and its stride size (consecutive elements / entry) for
	// both the data itself and the stride themselves.
	// Note that only INFO/FORMAT/FILTER fields have IDX fields. The
	// other fields do not require dictionary lookup to ascertain
	// their identity as they are guaranteed to be invariant.
	header_type* offsets;
	header_type* info_offsets;
	header_type* format_offsets;
	header_type* filter_offsets;

	// Bit-vectors of INFO/FORMAT/FILTER vectors of local IDX
	// patterns. These bit-vectors are used to quickly check
	// for the set membership of a given global and/or local IDX.
	// The bit-vectors internally holds the actual vector of IDX
	// for internal use. Construction of these bit-vectors are not
	// critical for basic functionality but critical for the
	// restoration of a bit-exact output sequence of fields.
	uint32_t n_info_patterns_allocated;
	uint32_t n_format_patterns_allocated;
	uint32_t n_filter_patterns_allocated;
	yon_blk_bv_pair* info_patterns;
	yon_blk_bv_pair* format_patterns;
	yon_blk_bv_pair* filter_patterns;

	// Supportive hash tables to permit the map from global
	// IDX fields to local IDX fields.
	map_type* info_map;
	map_type* format_map;
	map_type* filter_map;
	map_pattern_type* info_pattern_map;
	map_pattern_type* format_pattern_map;
	map_pattern_type* filter_pattern_map;
};

struct yon_blk_load_settings {
public:
	yon_blk_load_settings() : loaded_genotypes(false){}
	~yon_blk_load_settings() = default;

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

struct yon_vb_istats_obj {
	typedef yon_vb_istats_obj self_type;
	typedef yon1_dc_t  data_container_type;

	yon_vb_istats_obj(void) :
		cost_uncompressed(0),
		cost_strides(0),
		cost_compressed(0),
		cost_strides_compressed(0),
		cost_bcf(0)
	{}

	~yon_vb_istats_obj(){}

	self_type& operator+=(const self_type& other){
		this->cost_compressed += other.cost_compressed;
		this->cost_strides += other.cost_strides;
		this->cost_strides_compressed += other.cost_strides_compressed;
		this->cost_uncompressed += other.cost_uncompressed;
		this->cost_bcf += other.cost_bcf;
		return(*this);
	}

	self_type& operator+=(const data_container_type& container){
		this->cost_uncompressed += container.data_uncompressed.size();
		if(container.header.data_header.HasMixedStride())
			this->cost_strides += container.strides_uncompressed.size();

		this->cost_compressed   += container.data.size();
		if(container.header.data_header.HasMixedStride())
			this->cost_strides_compressed += container.strides.size();

		return(*this);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.cost_compressed << '\t' << entry.cost_uncompressed <<
				'\t' << entry.cost_strides_compressed << "\t" << entry.cost_strides <<
				'\t' << (entry.cost_compressed ? (double)(entry.cost_uncompressed)/(entry.cost_compressed) : 0) <<
				'\t' << entry.cost_bcf << '\t' << (entry.cost_compressed ? (double)entry.cost_bcf/(entry.cost_compressed) : 0);
		return(stream);
	}

	uint64_t cost_uncompressed;
	uint64_t cost_strides;
	uint64_t cost_compressed;
	uint64_t cost_strides_compressed;
	uint64_t cost_bcf;
};

class yon_vb_istats {
public:
	typedef yon_vb_istats self_type;
	typedef yon_vb_istats_obj    value_type;
    typedef std::size_t                   size_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;

public:
    yon_vb_istats() :
    	n_entries_(0),
		n_capacity_(0),
		entries_(nullptr)
    {
    }

    yon_vb_istats(const size_type start_capacity) :
    	n_entries_(0),
		n_capacity_(start_capacity),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {

    }

    ~yon_vb_istats(){
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~yon_vb_istats_obj();

		::operator delete[](static_cast<void*>(this->entries_));
    }

    yon_vb_istats(const self_type& other) :
    	n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {
    	for(uint32_t i = 0; i < this->size(); ++i)
    		new( &this->entries_[i] ) value_type( other.at(i) );
    }

    yon_vb_istats(self_type&& other) noexcept :
    	n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		entries_(nullptr)
    {
    	std::swap(this->entries_, other.entries_);
    }

    yon_vb_istats& operator=(self_type&& other) noexcept {
    	this->n_entries_ = other.n_entries_;
    	this->n_capacity_ = other.n_capacity_;
    	// Clear local entries if any.
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~yon_vb_istats_obj();

		::operator delete[](static_cast<void*>(this->entries_));
		this->entries_ = nullptr;

    	std::swap(this->entries_, other.entries_);
    	return(*this);
	}

    yon_vb_istats& operator=(const self_type& other) {
		this->n_entries_ = other.n_entries_;
		this->n_capacity_ = other.n_capacity_;
		// Clear local entries if any.
		for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~yon_vb_istats_obj();

		::operator delete[](static_cast<void*>(this->entries_));
		this->entries_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));

		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( other.at(i) );

		return(*this);
	}

	// Element access
	inline reference at(const size_type& position){ return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
	inline reference operator[](const size_type& position){ return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
	inline pointer data(void){ return(this->entries_); }
	inline const_pointer data(void) const{ return(this->entries_); }
	inline reference front(void){ return(this->entries_[0]); }
	inline const_reference front(void) const{ return(this->entries_[0]); }
	inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	self_type& operator+=(const self_type& other){
		while(other.size() > this->size()) this->resize();

		for(uint32_t i = 0; i < other.size(); ++i)
			this->at(i) += other[i];

		return(*this);
	}

	//
	inline self_type& operator+=(const const_reference entry){
		if(this->size() + 1 == this->n_capacity_)
			this->resize();


		this->entries_[this->n_entries_++] = entry;
		return(*this);
	}
	inline self_type& add(const const_reference entry){ return(*this += entry); }

	void resize(const uint32_t new_size){
		if(new_size == this->capacity()) return;
		if(new_size < this->capacity()){
			if(new_size < this->size()){
				for(uint32_t i = new_size; i < this->size(); ++i)
					(this->entries_ + i)->~yon_vb_istats_obj();
				this->n_entries_ = new_size;
			}
			return;
		}

		pointer temp = this->entries_;

		this->n_capacity_ = new_size;
		this->entries_    = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));

		// Lift over values from old addresses
		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( temp[i] );

		if(temp != nullptr){
			for(std::size_t i = 0; i < this->size(); ++i)
				(temp + i)->~yon_vb_istats_obj();

			::operator delete[](static_cast<void*>(temp));
		}
	}

	inline void resize(void){ return(this->resize(this->capacity()*2)); }

	void Allocate(const uint32_t entries){
		// Delete previous entries.
		if(this->entries_ != nullptr){
			for(std::size_t i = 0; i < this->size(); ++i)
				(this->entries_ + i)->~yon_vb_istats_obj();
			::operator delete[](static_cast<void*>(this->entries_));

			this->entries_ = nullptr;
		}

		if(entries > this->capacity()) this->resize(entries);

		this->n_entries_ = entries;
		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( );
	}

private:
	size_type n_entries_;
	size_type n_capacity_;
	pointer   entries_;
};

/**
 * Primary Tachyon block object: stores containers of data and
 * provides encapsulated and abstracted access to its
 * contents.
 */
class yon1_vb_t {
public:
	typedef yon1_vb_t                    self_type;
	typedef yon1_dc_t                   container_type;
	typedef yon_vb_hdr              block_header_type;
	typedef yon_vb_ftr              block_footer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef yon_vb_istats import_stats_type;
	typedef yon_dc_hdr             offset_type;
	typedef DataBlockSettings               block_settings_type;

public:
	yon1_vb_t();
	yon1_vb_t(const uint16_t n_info, const uint16_t n_format);
	~yon1_vb_t();
	yon1_vb_t(const self_type& other);
	yon1_vb_t(self_type&& other) noexcept;
	yon1_vb_t& operator=(const self_type& other);
	yon1_vb_t& operator=(self_type&& other) noexcept;

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
	 * @return             Returns TRUE upon success or FALSE otherwise
	 */
	bool write(std::ostream& stream);

	/**<
	 * Add the data from a MetaEntry object to this block. Internally
	 * performs all the operations required to transfer each MetaEntry
	 * field into the correct destination with the correct encodings.
	 * This is the preferred way of storing a MetaEntry.
	 * @param meta_entry Input MetaEntry object to be stored.
	 * @return           Returns TRUE upon success or FALSE otherwise.
	 */
	bool operator+=(yon1_vnt_t& rcd);

	/**<
	 * Wrapper to add a new DataContainer to this yon1_vb_t container.
	 * @param dc
	 * @param dst_containers
	 * @param AddWrapper
	 * @param StreamFieldLookup
	 * @return
	 */
	bool AddWrapper(yon1_dc_t& dc,
					const int32_t idx,
	                container_type* dst_containers,
					uint32_t(yon_vb_ftr::*AddWrapper)(const uint32_t),
	                int32_t(yon1_vb_t::*StreamFieldLookup)(const uint32_t) const)
	{
		if(idx < 0){
			std::cerr << "illegal idx" << std::endl;
			return false;
		}
		dc.header.data_header.global_key = idx;

		// 1) Search for global idx
		int field_exists = (this->*StreamFieldLookup)(dc.GetIdx());
		if(field_exists >= 0){
			//assert(dst_containers[field_exists].GetDataPrimitiveType() == dc.GetDataPrimitiveType());
			dst_containers[field_exists] += dc;
			//assert(dst_containers[field_exists].header.data_header.stride == dc.header.data_header.stride);
			//std::cerr << "merging: " << dst_containers[field_exists].header.data_header.stride << " and " << dc.header.data_header.stride << std::endl;
			//std::cerr << "field exists at " << field_exists << " for id " << (uint32_t)dc.GetIdx() << ": adding now " << dst_containers[field_exists].GetSizeUncompressed() << std::endl;
		} else {
			uint32_t offset = (this->footer.*AddWrapper)(dc.GetIdx());
			//std::cerr << "field does not exist. adding. " << offset  << std::endl;
			dst_containers[offset] = dc;
		}

		return true;
	}

	inline bool AddInfo(yon1_dc_t& dc, const YonInfo* info){
		return(this->AddWrapper(dc,
		       info->idx,
		       this->info_containers,
		       &block_footer_type::AddInfo,
		       &self_type::GetInfoPosition));
	}

	inline bool AddFormat(yon1_dc_t& dc, const YonFormat* fmt){
		return(this->AddWrapper(dc,
		                        fmt->idx,
		                        this->format_containers,
		                        &block_footer_type::AddFormat,
		                        &self_type::GetFormatPosition));
	}

	self_type& AddMore(yon1_vnt_t& rec){
		if(this->header.n_variants == 0){
			this->header.contig_id    = rec.rid;
			this->header.min_position = rec.pos;
			this->header.max_position = rec.pos;
		}

		++this->header.n_variants;
		this->header.max_position = rec.pos;

		// TODO: MAJOR WARNING: NO GUARANTEE THAT DATA ADDED TO THE BYTE STREAMS
		//       HAVE THE SAME PRIMITIVE TYPE AS IS __REQUIRED__ IN TACHYON.
		std::vector<int> id_list(rec.n_info);
		if(rec.n_info){
			for(int i = 0; i < rec.n_info; ++i){
				id_list[i] = rec.info_hdr[i]->idx;
				// 1 add info, format, and filter
				yon1_dc_t dc(std::move(rec.info[i]->ToDataContainer()));

				if(rec.info_hdr[i]->yon_type == YON_VCF_HEADER_INTEGER){
					dc.header.data_header.controller.type = YON_TYPE_32B;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.info_hdr[i]->yon_type == YON_VCF_HEADER_FLOAT){
					dc.header.data_header.controller.type = YON_TYPE_FLOAT;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.info_hdr[i]->yon_type == YON_VCF_HEADER_STRING || rec.info_hdr[i]->yon_type == YON_VCF_HEADER_CHARACTER){
					dc.header.data_header.controller.type = YON_TYPE_CHAR;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.info_hdr[i]->yon_type == YON_VCF_HEADER_FLAG){
					dc.header.data_header.controller.type = YON_TYPE_BOOLEAN;
					dc.header.data_header.controller.signedness = false;
				} else {
					std::cerr << "unknwon type" << std::endl;
					exit(1);
				}

				if(this->AddInfo(dc, rec.info_hdr[i]) == false){
					std::cerr << "failed to add info: " << rec.info_hdr[i]->id << std::endl;
				}
			}

			rec.info_pid = this->AddInfoPattern(id_list);
			//std::cerr << "info_pid=" << rec.info_pid << std::endl;
		} else rec.info_pid = -1;

		if(rec.n_fmt){
			id_list.resize(rec.n_fmt);
			for(int i = 0; i < rec.n_fmt; ++i){
				id_list[i] = rec.fmt_hdr[i]->idx;
				yon1_dc_t dc(std::move(rec.fmt[i]->ToDataContainer()));

				if(rec.fmt_hdr[i]->yon_type == YON_VCF_HEADER_INTEGER){
					dc.header.data_header.controller.type = YON_TYPE_32B;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.fmt_hdr[i]->yon_type == YON_VCF_HEADER_FLOAT){
					dc.header.data_header.controller.type = YON_TYPE_FLOAT;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.fmt_hdr[i]->yon_type == YON_VCF_HEADER_STRING || rec.fmt_hdr[i]->yon_type == YON_VCF_HEADER_CHARACTER){
					dc.header.data_header.controller.type = YON_TYPE_CHAR;
					dc.header.data_header.controller.signedness = true;
				} else if(rec.fmt_hdr[i]->yon_type == YON_VCF_HEADER_FLAG){
					dc.header.data_header.controller.type = YON_TYPE_BOOLEAN;
					dc.header.data_header.controller.signedness = false;
				} else {
					std::cerr << "unknwon type" << std::endl;
					exit(1);
				}

				if(this->AddFormat(dc, rec.fmt_hdr[i]) == false){
					std::cerr << "failed to add format: " << rec.fmt_hdr[i]->id << std::endl;
				}
			}

			rec.fmt_pid = this->AddFormatPattern(id_list);
			//std::cerr << "fmt_pid=" << rec.fmt_pid << std::endl;
		} else rec.fmt_pid = -1;

		if(rec.n_flt){
			id_list.resize(rec.n_flt);
			for(int i = 0; i < rec.n_flt; ++i){
				id_list[i] = rec.flt_hdr[i]->idx;
				this->AddFilter(id_list[i]);
			}
			rec.flt_pid = this->AddFilterPattern(id_list);
		} else rec.fmt_pid = -1;

		return(*this);
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
	void UpdateContainers(const uint32_t n_sampless);

	/**<
	 * Determine compressed block-size. Execute this function prior to writing a
	 * block.
	 * @return Returns the total number of bytes that will be written.
	 */
	uint64_t GetCompressedSize(void) const;
	uint64_t GetUncompressedSize(void) const;


	inline void PackFooter(void){
		this->footer_support.reset();
		this->footer_support.data_uncompressed << this->footer;
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
		if(this->header.controller.has_gt && this->header.controller.has_gt_permuted){
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

	/**<
	 *
	 * @param stats_basic
	 * @param stats_info
	 * @param stats_format
	 */
	void UpdateOutputStatistics(import_stats_type& stats_basic,
								import_stats_type& stats_info,
								import_stats_type& stats_format);

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
	 * Move over pair of headers from a data container to a block footer
	 * @param offset    Destination header in footer
	 * @param container Target container hosting the header
	 */
	inline static void UpdateHeader(offset_type& offset, const container_type& container){
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
	static void UpdateHeader(offset_type& offset,
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
	 * Write a src container to a dst output stream.
	 * @param stream         Dst output stream.
	 * @param offset         Dst offset header in the footer to update.
	 * @param container      Src container to write.
	 * @param virtual_offset Virtual file/stream offset at the start of this container.
	 */
	void static WriteContainer(std::ostream& stream,
	                           offset_type& offset,
	                           const container_type& container,
	                           const uint32_t virtual_offset)
	{
		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
			return(yon1_vb_t::WriteContainerEncrypted(stream, offset, container, virtual_offset));

		assert(offset.data_header.offset == virtual_offset);
		yon1_vb_t::UpdateHeader(offset, container, virtual_offset);
		assert(container.data.size() == offset.data_header.cLength);
		stream << container;
	}

	/**<
	 * Write a src container to a dst output stream when the src container
	 * is partially/completely encrypted.
	 * @param stream
	 * @param offset
	 * @param container
	 * @param virtual_offset
	 */
	static void WriteContainerEncrypted(std::ostream& stream,
	                                    offset_type& offset,
	                                    const container_type& container,
	                                    const uint32_t virtual_offset)
	{
		yon1_vb_t::UpdateHeader(offset, container, virtual_offset);
		assert(container.data.size() == offset.data_header.eLength);
		// Encrypted data is concatenated: write only data buffer
		stream.write(container.data.data(), container.data.size());
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

//
class yon_cont_ref_iface {
public:
	yon_cont_ref_iface() : is_uniform_(false), n_offset_(0), n_elements_(0), b_data_(0), data_(nullptr){}
	yon_cont_ref_iface(char* data, uint64_t l_data) : is_uniform_(false), n_offset_(0), n_elements_(0), b_data_(l_data), data_(data){}
	virtual ~yon_cont_ref_iface(){
		// do not delete data. it is not owned by this
	}

	// next()
	// -> stride of response + response pointer
	virtual bool next() =0;
	virtual int32_t GetInt32(const uint32_t p) =0;

public:
	bool is_uniform_;
	uint32_t n_offset_;
	uint32_t n_elements_;
	uint64_t b_data_;
	char* data_;
};

template <class return_ptype>
class yon_cont_ref : public yon_cont_ref_iface {
public:
	yon_cont_ref() : variants_(nullptr){}
	yon_cont_ref(char* data, uint64_t l_data) :
		yon_cont_ref_iface(data, l_data),
		variants_(reinterpret_cast<return_ptype*>(this->data_))
	{
		this->n_elements_ = l_data / sizeof(return_ptype);
		assert(l_data % sizeof(return_ptype) == 0);
	}
	~yon_cont_ref(){}

	bool next(){
		if(this->is_uniform_) return true;
		if(this->n_offset_ == this->n_elements_) return false;
		++this->n_offset_;
		return true;
	}

	inline return_ptype& operator[](const uint32_t p){ return(this->variants_[p]); }
	inline const return_ptype& operator[](const uint32_t p) const{ return(this->variants_[p]); }
	inline int32_t GetInt32(const uint32_t p){ return((int32_t)this->variants_[p]); }

public:
	return_ptype* variants_;
};

class yon1_vc_t {
public:
	typedef yon1_dc_t dc_type;
	typedef yon1_vc_t self_type;
    typedef std::size_t       size_type;
    typedef yon1_vnt_t        value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
	yon1_vc_t() : n_variants_(0), n_capacity_(0), variants_(nullptr){}
	yon1_vc_t(const uint32_t size) :
		n_variants_(0), n_capacity_(size),
		//variants_(new yon1_vnt_t[size])
		variants_(static_cast<yon1_vnt_t*>(::operator new[](this->n_capacity_*sizeof(yon1_vnt_t))))
	{

	}

	yon1_vc_t(yon1_vb_t& variant_block, const VariantHeader& header) :
		n_variants_(variant_block.header.n_variants), n_capacity_(variant_block.header.n_variants + 64),
		//variants_(new yon1_vnt_t[size])
		variants_(static_cast<yon1_vnt_t*>(::operator new[](this->n_capacity_*sizeof(yon1_vnt_t))))
	{
		if(this->n_variants_){
			for(int i = 0; i < this->n_variants_; ++i)
				new( &this->variants_[i] ) yon1_vnt_t( );

			this->Build(variant_block, header);
		}
	}

	~yon1_vc_t(){
		//delete [] variants_;
		for(std::size_t i = 0; i < this->size(); ++i)
			((this->variants_ + i)->~yon1_vnt_t)();

		::operator delete[](static_cast<void*>(this->variants_));
	}

	yon1_vc_t& operator+=(const yon1_vnt_t& rec){
		if(this->variants_ == nullptr || this->capacity() == 0)
			this->reserve();

		if(this->size() == this->capacity())
			this->reserve();

		new( &this->variants_[this->n_variants_] ) yon1_vnt_t( rec );
		++this->n_variants_;

		return(*this);
	}

	const size_t size(void) const { return(this->n_variants_); }
	const size_t capacity(void) const { return(this->n_capacity_); }

	void clear(void){
		for(std::size_t i = 0; i < this->size(); ++i)
			((this->variants_ + i)->~yon1_vnt_t)();

		this->n_variants_ = 0;
		this->block_.clear();
	}

	inline void reserve(void){ this->reserve(n_capacity_ + 1000); }
	void reserve(const uint32_t new_size);
	void resize(const uint32_t new_size);

	bool Build(yon1_vb_t& variant_block, const VariantHeader& header);

	/**< @brief Reads a target Tachyon block from disk.
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity in this function.
	 * @param stream    Src stream handle.
	 * @param header    Src global header.
	 * @param settings  Settings record describing reading parameters.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	inline bool ReadBlock(std::ifstream& stream, const VariantHeader& header, DataBlockSettings& settings){
		return(this->block_.read(stream, settings, header));
	}

	bool PrepareWritableBlock(const uint32_t n_samples);

	index::VariantIndexEntry GetIndexEntry(void){
		index::VariantIndexEntry index_entry;
		//index_entry.blockID     = block_id;
		index_entry.contig_id    = this->block_.header.contig_id;
		index_entry.min_position = this->block_.header.min_position;
		index_entry.max_position = this->block_.header.max_position;
		index_entry.n_variants   = this->block_.header.n_variants;
		return(index_entry);
	}

	 // Element access
	inline reference at(const size_type& position){ return(this->variants_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->variants_[position]); }
	inline reference operator[](const size_type& position){ return(this->variants_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->variants_[position]); }
	inline pointer data(void){ return(this->variants_); }
	inline const_pointer data(void) const{ return(this->variants_); }
	inline reference front(void){ return(this->variants_[0]); }
	inline const_reference front(void) const{ return(this->variants_[0]); }
	inline reference back(void){ return(this->variants_[(n_variants_ == 0 ? 0 : n_variants_ - 1)]); }
	inline const_reference back(void) const{ return(this->variants_[(n_variants_ == 0 ? 0 : n_variants_ - 1)]); }

	// Iterator
	inline iterator begin(){ return iterator(&this->variants_[0]); }
	inline iterator end(){ return iterator(&this->variants_[this->n_variants_]); }
	inline const_iterator begin() const{ return const_iterator(&this->variants_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->variants_[this->n_variants_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->variants_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->variants_[this->n_variants_]); }

private:
	bool PermuteOrder(const yon1_vb_t& variant_block);

	bool AddInfo(yon1_vb_t& variant_block, const VariantHeader& header);
	bool AddFilter(yon1_vb_t& variant_block, const VariantHeader& header);
	bool AddFormat(yon1_vb_t& variant_block, const VariantHeader& header);
	bool AddInfoWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool AddFormatWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool InfoSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool InfoSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool AddGenotypes(yon1_vb_t& block, const VariantHeader& header);

	template <class T>
	bool AddBaseInteger(dc_type& container, void(yon1_vnt_t::*fnc)(const T v)){
		if(container.data_uncompressed.size() == 0)
			return false;

		assert(container.header.HasMixedStride() == false);

		yon_cont_ref_iface* it = nullptr;
		switch(container.header.data_header.controller.type){
		case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
		case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
		case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
		case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
		}
		assert(it != nullptr);

		// If data is uniform.
		if(container.header.data_header.controller.uniform){
			it->is_uniform_ = true;

			for(int i = 0; i < this->n_variants_; ++i){
				(variants_[i].*fnc)(it->GetInt32(0));
			}
		} else {
			assert(this->n_variants_ == it->n_elements_);
			for(int i = 0; i < this->n_variants_; ++i){
				(variants_[i].*fnc)(it->GetInt32(i));
			}
		}
		delete it;
		return true;
	}

	bool FormatSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool FormatSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride);

	inline bool AddContigs(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetChromosome)); }
	inline bool AddController(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetController)); }
	inline bool AddPositions(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetPosition)); }

	bool AddQuality(dc_type& container);
	bool AddRefAlt(dc_type& container);

	inline bool AddFilterIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFilterPId)); }
	inline bool AddFormatIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFormatPId)); }
	inline bool AddInfoIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetInfoPId)); }
	inline bool AddPloidy(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetBasePloidy)); }

	bool AddAlleles(dc_type& container);
	bool AddNames(dc_type& container);

public:
	yon1_vb_t block_;
	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	//yon_gt_rcd** gt_exp;

	size_t n_variants_;
	size_t n_capacity_;
	yon1_vnt_t* variants_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::InfoSetup(dc_type& container,
                                 const VariantHeader& header,
                                 const std::vector<bool>& matches)
{
	if(container.strides_uncompressed.size() == 0)
		return false;

	yon_cont_ref_iface* it = nullptr;
	switch(container.header.stride_header.controller.type){
	case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	}
	assert(it != nullptr);

	uint32_t current_offset = 0;
	uint32_t stride_offset = 0;

	for(int i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].info_pid == -1){
			continue;
		} else if(matches[this->variants_[i].info_pid]){
			this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, it->GetInt32(stride_offset));
			this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
			current_offset += it->GetInt32(stride_offset) * sizeof(intrinsic_ptype);
			++stride_offset;
		}
	}
	assert(current_offset == container.data_uncompressed.size());
	delete it;
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::InfoSetup(dc_type& container,
                                 const VariantHeader& header,
                                 const std::vector<bool>& matches,
                                 const uint32_t stride_size)
{
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, 0, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype);
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::FormatSetup(dc_type& container,
                                   const VariantHeader& header,
                                   const std::vector<bool>& matches)
{
	if(container.strides_uncompressed.size() == 0)
		return false;

	yon_cont_ref_iface* it = nullptr;
	switch(container.header.stride_header.controller.type){
	case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	}
	assert(it != nullptr);

	uint32_t current_offset = 0;
	uint32_t stride_offset = 0;

	for(int i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].fmt_pid == -1){
			continue;
		} else if(matches[this->variants_[i].fmt_pid]){
			this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), it->GetInt32(stride_offset));
			this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
			current_offset += it->GetInt32(stride_offset) * sizeof(intrinsic_ptype) * header.GetNumberSamples();
			++stride_offset;
		}
	}
	assert(current_offset == container.data_uncompressed.size());
	delete it;
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::FormatSetup(dc_type& container,
                                   const VariantHeader& header,
                                   const std::vector<bool>& matches,
                                   const uint32_t stride_size)
{
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, 0, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype) * header.GetNumberSamples();
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

}



#endif /* CONTAINERS_VARIANT_CONTAINERS_H_ */
