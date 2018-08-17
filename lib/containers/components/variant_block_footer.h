#ifndef CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_
#define CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_

#include <unordered_map>

#include "data_container_header.h"
#include "io/basic_buffer.h"

#include "third_party/xxhash/xxhash.h"

namespace tachyon {

#define YON_BLK_N_STATIC   25// Total number of invariant headers
#define YON_BLK_PPA         0 // Sample permutation array
#define YON_BLK_CONTIG      1
#define YON_BLK_POSITION    2
#define YON_BLK_REFALT      3
#define YON_BLK_CONTROLLER  4 // Set memberships
#define YON_BLK_QUALITY     5
#define YON_BLK_NAMES       6
#define YON_BLK_ALLELES     7
#define YON_BLK_ID_INFO     8
#define YON_BLK_ID_FORMAT   9
#define YON_BLK_ID_FILTER  10
#define YON_BLK_GT_INT8    11 // Run-length encoded genotypes
#define YON_BLK_GT_INT16   12
#define YON_BLK_GT_INT32   13
#define YON_BLK_GT_INT64   14
#define YON_BLK_GT_S_INT8  15 // Standard encoded genotypes
#define YON_BLK_GT_S_INT16 16
#define YON_BLK_GT_S_INT32 17
#define YON_BLK_GT_S_INT64 18
#define YON_BLK_GT_N_INT8  19 // Standard encoded genotypes
#define YON_BLK_GT_N_INT16 20
#define YON_BLK_GT_N_INT32 21
#define YON_BLK_GT_N_INT64 22
#define YON_BLK_GT_SUPPORT 23 // Genotype support
#define YON_BLK_GT_PLOIDY  24 // Genotype ploidy

#define YON_BLK_BV_PPA         1 << YON_BLK_PPA
#define YON_BLK_BV_CONTIG      1 << YON_BLK_CONTIG
#define YON_BLK_BV_POSITION    1 << YON_BLK_POSITION
#define YON_BLK_BV_REFALT      1 << YON_BLK_REFALT
#define YON_BLK_BV_CONTROLLER  1 << YON_BLK_CONTROLLER
#define YON_BLK_BV_QUALITY     1 << YON_BLK_QUALITY
#define YON_BLK_BV_NAMES       1 << YON_BLK_NAMES
#define YON_BLK_BV_ALLELES     1 << YON_BLK_ALLELES
#define YON_BLK_BV_ID_INFO     1 << YON_BLK_ID_INFO
#define YON_BLK_BV_ID_FORMAT   1 << YON_BLK_ID_FORMAT
#define YON_BLK_BV_ID_FILTER   1 << YON_BLK_ID_FILTER
#define YON_BLK_BV_GT_INT8     1 << YON_BLK_GT_INT8
#define YON_BLK_BV_GT_INT16    1 << YON_BLK_GT_INT16
#define YON_BLK_BV_GT_INT32    1 << YON_BLK_GT_INT32
#define YON_BLK_BV_GT_INT64    1 << YON_BLK_GT_INT64
#define YON_BLK_BV_GT_S_INT8   1 << YON_BLK_GT_S_INT8
#define YON_BLK_BV_GT_S_INT16  1 << YON_BLK_GT_S_INT16
#define YON_BLK_BV_GT_S_INT32  1 << YON_BLK_GT_S_INT32
#define YON_BLK_BV_GT_S_INT64  1 << YON_BLK_GT_S_INT64
#define YON_BLK_BV_GT_N_INT8   1 << YON_BLK_GT_N_INT8
#define YON_BLK_BV_GT_N_INT16  1 << YON_BLK_GT_N_INT16
#define YON_BLK_BV_GT_N_INT32  1 << YON_BLK_GT_N_INT32
#define YON_BLK_BV_GT_N_INT64  1 << YON_BLK_GT_N_INT64
#define YON_BLK_BV_GT_SUPPORT  1 << YON_BLK_GT_SUPPORT
#define YON_BLK_BV_GT_PLOIDY   1 << YON_BLK_GT_PLOIDY

#define YON_BLK_BV_INFO        1 << (YON_BLK_N_STATIC)
#define YON_BLK_BV_FORMAT      1 << (YON_BLK_N_STATIC + 1)
#define YON_BLK_BV_GT          ((YON_BLK_BV_GT_INT8)|(YON_BLK_BV_GT_INT16)|(YON_BLK_BV_GT_INT32)|(YON_BLK_BV_GT_INT64)|(YON_BLK_BV_GT_S_INT8)|(YON_BLK_BV_GT_S_INT16)|(YON_BLK_BV_GT_S_INT32)|(YON_BLK_BV_GT_S_INT64)|(YON_BLK_BV_GT_N_INT8)|(YON_BLK_BV_GT_N_INT16)|(YON_BLK_BV_GT_N_INT32)|(YON_BLK_BV_GT_N_INT64)|(YON_BLK_BV_GT_SUPPORT)|(YON_BLK_BV_GT_PLOIDY))

namespace containers {

struct yon_blk_bv_pair {
	yon_blk_bv_pair() : l_bytes(0), bit_bytes(nullptr){}
	~yon_blk_bv_pair(){ delete [] this->bit_bytes; }

	void clear(void){
		this->pattern.clear();
		this->l_bytes = 0;
		delete [] this->bit_bytes;
		this->bit_bytes = nullptr;
	}

	// Bit access
	inline bool operator[](const U32 position) const{ return((this->bit_bytes[position / 8] & (1 << (position % 8))) >> (position % 8)); }


	// Given the total number of fields allocate ceil(n_total_fields/8)
	// bytes for the base array.
	void Build(const U32 n_total_fields, const std::unordered_map<U32, U32>* local_map){
		if(this->pattern.size() == 0) return;
		assert(local_map != nullptr);

		// Determine the required width in bytes of the bit-vector
		BYTE bitvector_width = ceil((float)(n_total_fields+1)/8);

		// Allocate new bit-vectors
		delete [] this->bit_bytes;
		this->l_bytes = bitvector_width;
		this->bit_bytes = new uint8_t[bitvector_width];

		// Cycle over pattern size
		for(U32 i = 0; i < this->pattern.size(); ++i){
			std::unordered_map<U32, U32>::const_iterator it = local_map->find(this->pattern[i]);
			assert(it != local_map->end());

			// Map from absolute key to local key.
			U32 local_key = it->second;
			assert(local_key <= n_total_fields);

			// Set bit at local key position
			this->bit_bytes[local_key/8] |= 1 << (local_key % 8);
		}
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const yon_blk_bv_pair& entry){
		io::SerializePrimitive(entry.l_bytes, buffer);
		buffer += (U32)entry.pattern.size();
		for(U32 i = 0; i < entry.pattern.size(); ++i)
			io::SerializePrimitive(entry.pattern[i], buffer);

		for(U32 i = 0; i < entry.l_bytes; ++i)
			io::SerializePrimitive(entry.bit_bytes[i], buffer);


		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, yon_blk_bv_pair& entry){
		entry.pattern.clear();
		io::DeserializePrimitive(entry.l_bytes, buffer);
		U32 l_vector;
		buffer >> l_vector;
		//entry.pattern.resize(l_vector);
		for(U32 i = 0; i < l_vector; ++i){
			int temp;
			io::DeserializePrimitive(temp, buffer);
			//entry.pattern[i] = temp;
			entry.pattern.push_back(temp);
		}

		entry.bit_bytes = new BYTE[entry.l_bytes];
		for(U32 i = 0; i < entry.l_bytes; ++i)
			io::DeserializePrimitive(entry.bit_bytes[i], buffer);

		return(buffer);
	}

	yon_blk_bv_pair& operator=(const yon_blk_bv_pair& other){
		delete [] this->bit_bytes;
		this->pattern = other.pattern;
		this->l_bytes = other.l_bytes;
		this->bit_bytes = new uint8_t[this->l_bytes];
		memcpy(this->bit_bytes, other.bit_bytes, this->l_bytes);
		return(*this);
	}

	yon_blk_bv_pair& operator=(yon_blk_bv_pair&& other) noexcept{
		if (this == &other){
			// take precautions against self-moves
			return *this;
		}

		delete [] this->bit_bytes; this->bit_bytes = nullptr;
		std::swap(this->bit_bytes, other.bit_bytes);
		this->pattern = std::move(other.pattern);
		this->l_bytes = other.l_bytes;
		return(*this);
	}

public:
	std::vector<int> pattern;
	uint8_t l_bytes;
	uint8_t* bit_bytes;
};

// It is possible of getting mapping local indices to global IDX
// for either  FILTER/FORMAT/INFO fields by iterating over the
// relevant DataContainerHeader structures and tracking the incremental
// position they occur at (local IDX) and read the global IDX in the
// header structure.
struct VariantBlockFooter {
public:
	typedef VariantBlockFooter        self_type;
	typedef DataContainerHeader       header_type;
	typedef std::unordered_map<U32, U32>    map_type;
	typedef std::unordered_map<U64, U32>    map_pattern_type;

public:
	VariantBlockFooter();
	~VariantBlockFooter();
	VariantBlockFooter(const self_type& other);
	VariantBlockFooter(self_type&& other) noexcept;
	VariantBlockFooter& operator=(const self_type& other);
	VariantBlockFooter& operator=(self_type&& other) noexcept;

	void reset(void);
	void resetTables(void);

	// Allocate offset vectors
	void AllocateInfoHeaders(const U32 n_info_streams);
	void AllocateFormatHeaders(const U32 n_format_streams);
	void AllocateFilterHeaders(const U32 n_filter_streams);

	/**<
	 * Wrapper function for allocating new offset objects for INFO,
	 * FORMAT, and FILTER patterns and streams
	 * @param n_info_streams   Number of unique info streams
	 * @param n_format_streams Number of unique format streams
	 * @param n_filter_streams Number of unique filter streams
	 */
	void AllocateHeaders(const U32 n_info_streams,
		                        const U32 n_format_streams,
		                        const U32 n_filter_streams);
	bool ConstructInfoBitVector(std::unordered_map<U32,U32>* pattern_map);
	bool ConstructFormatBitVector(std::unordered_map<U32,U32>* pattern_map);
	bool ConstructFilterBitVector(std::unordered_map<U32,U32>* pattern_map);

	U32 AddPatternWrapper(const std::vector<int>& pattern,
	                      map_pattern_type* pattern_map,
	                      yon_blk_bv_pair* bv_pairs,
	                      U16& stream_counter);
	U32 AddInfoPattern(const std::vector<int>& pattern);
	U32 AddFormatPattern(const std::vector<int>& pattern);
	U32 AddFilterPattern(const std::vector<int>& pattern);

	// This wrapper adds patterns to the hash map when the data has
	// already been loaded. This occurs when loading an object from
	// disk/buffer.
	U32 UpdatePatternWrapper(const std::vector<int>& pattern,
	                         map_pattern_type* pattern_map,
	                         const U16& stream_counter);

	U32 UpdateInfoPattern(const std::vector<int>& pattern, const U16 pattern_id);
	U32 UpdateFormatPattern(const std::vector<int>& pattern, const U16 pattern_id);
	U32 UpdateFilterPattern(const std::vector<int>& pattern, const U16 pattern_id);

	void Finalize(void);
	bool BuildMaps(void);
	bool BuildPatternMaps(void);

	U32 UpdateOffsetMapWrapper(const header_type& offset, map_type* map, const U16& stream_counter);
	U32 UpdateInfo(const header_type& offset, const U16 position);
	U32 UpdateFormat(const header_type& offset, const U16 position);
	U32 UpdateFilter(const header_type& offset, const U16 position);

	U32 AddStreamWrapper(const U32 id, map_type* map, header_type*& offsets, U16& stream_counter);
	U32 AddInfo(const U32 id);
	U32 AddFormat(const U32 id);
	U32 AddFilter(const U32 id);

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

		XXH_errorcode const resetResult = XXH64_reset(state, 71236251);
		if (resetResult == XXH_ERROR) abort();

		for(U32 i = 0; i < id_vector.size(); ++i){
			XXH_errorcode const addResult = XXH64_update(state, (const void*)&id_vector[i], sizeof(int));
			if (addResult == XXH_ERROR) abort();
		}

		U64 hash = XXH64_digest(state);
		XXH64_freeState(state);

		return hash;
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& entry);

public:
	// Utility members. l_*_bitvector stores the byte-length
	// of each of the bit-vectors described below.
	// Not written or read from disk. Used internally only
	BYTE l_info_bitvector;
	BYTE l_format_bitvector;
	BYTE l_filter_bitvector;

	// Critical values used to track the number of data streams
	// that is available for each possible type. The n_*_streams
	// fields corresponds to the number of different byte streams
	// that are set in this block. The n_*_patterns corresponds
	// to the number of unique vectors of field identifiers that
	// occurred in the block.
	U16 n_info_streams; // streams
	U16 n_format_streams;
	U16 n_filter_streams;
	U16 n_info_patterns; // patterns
	U16 n_format_patterns;
	U16 n_filter_patterns;

	// Header structures corresponds critical information regarding
	// the global IDX and virtual byte offset to the start of each
	// compressed and possibly encrypted byte stream. In addition,
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
	U32 n_info_patterns_allocated;
	U32 n_format_patterns_allocated;
	U32 n_filter_patterns_allocated;
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

}
}



#endif /* CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_ */
