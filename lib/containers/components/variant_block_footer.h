#ifndef CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_
#define CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_

#include <unordered_map>

#include "data_container_header.h"
#include "io/basic_buffer.h"

#include "third_party/xxhash/xxhash.h"

namespace tachyon {

const std::vector<std::string> YON_BLK_PRINT_NAMES = {
		"PPA","CONTIG","POSITION","REFALT","CONTROLLER","QUALITY","NAMES",
		"ALLELES","ID_INFO","ID_FORMAT","ID_FILTER",
		"GT_INT8","GT_INT16","GT_INT32","GT_INT64",
		"GT_S_INT8","GT_S_INT16","GT_S_INT32","GT_S_INT64",
		"GT_N_INT8","GT_N_INT16","GT_N_INT32","GT_N_INT64",
		"GT_SUPPORT","GT_PLOIDY"};

/**<
 * These definitions correspond to the array offsets for the
 * invariant containers in the VariantBlock/VariantBlockFooter.
 * These correspond to the core components for site information
 * and internals such as CONTROLLER, ID_*, and GT_SUPPORT and
 * GT_PLOIDY.
 */
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

#define YON_BLK_BV_PPA         (1 << YON_BLK_PPA)
#define YON_BLK_BV_CONTIG      (1 << YON_BLK_CONTIG)
#define YON_BLK_BV_POSITION    (1 << YON_BLK_POSITION)
#define YON_BLK_BV_REFALT      (1 << YON_BLK_REFALT)
#define YON_BLK_BV_CONTROLLER  (1 << YON_BLK_CONTROLLER)
#define YON_BLK_BV_QUALITY     (1 << YON_BLK_QUALITY)
#define YON_BLK_BV_NAMES       (1 << YON_BLK_NAMES)
#define YON_BLK_BV_ALLELES     (1 << YON_BLK_ALLELES)
#define YON_BLK_BV_ID_INFO     (1 << YON_BLK_ID_INFO)
#define YON_BLK_BV_ID_FORMAT   (1 << YON_BLK_ID_FORMAT)
#define YON_BLK_BV_ID_FILTER   (1 << YON_BLK_ID_FILTER)
#define YON_BLK_BV_GT_INT8     (1 << YON_BLK_GT_INT8)
#define YON_BLK_BV_GT_INT16    (1 << YON_BLK_GT_INT16)
#define YON_BLK_BV_GT_INT32    (1 << YON_BLK_GT_INT32)
#define YON_BLK_BV_GT_INT64    (1 << YON_BLK_GT_INT64)
#define YON_BLK_BV_GT_S_INT8   (1 << YON_BLK_GT_S_INT8)
#define YON_BLK_BV_GT_S_INT16  (1 << YON_BLK_GT_S_INT16)
#define YON_BLK_BV_GT_S_INT32  (1 << YON_BLK_GT_S_INT32)
#define YON_BLK_BV_GT_S_INT64  (1 << YON_BLK_GT_S_INT64)
#define YON_BLK_BV_GT_N_INT8   (1 << YON_BLK_GT_N_INT8)
#define YON_BLK_BV_GT_N_INT16  (1 << YON_BLK_GT_N_INT16)
#define YON_BLK_BV_GT_N_INT32  (1 << YON_BLK_GT_N_INT32)
#define YON_BLK_BV_GT_N_INT64  (1 << YON_BLK_GT_N_INT64)
#define YON_BLK_BV_GT_SUPPORT  (1 << YON_BLK_GT_SUPPORT)
#define YON_BLK_BV_GT_PLOIDY   (1 << YON_BLK_GT_PLOIDY)

#define YON_BLK_BV_INFO        (1 << (YON_BLK_N_STATIC))
#define YON_BLK_BV_FORMAT      (1 << (YON_BLK_N_STATIC + 1))
#define YON_BLK_BV_GT          ((YON_BLK_BV_GT_INT8)|(YON_BLK_BV_GT_INT16)|(YON_BLK_BV_GT_INT32)|(YON_BLK_BV_GT_INT64)|(YON_BLK_BV_GT_S_INT8)|(YON_BLK_BV_GT_S_INT16)|(YON_BLK_BV_GT_S_INT32)|(YON_BLK_BV_GT_S_INT64)|(YON_BLK_BV_GT_N_INT8)|(YON_BLK_BV_GT_N_INT16)|(YON_BLK_BV_GT_N_INT32)|(YON_BLK_BV_GT_N_INT64)|(YON_BLK_BV_GT_SUPPORT)|(YON_BLK_BV_GT_PLOIDY))

struct yon_blk_bv_pair {
	yon_blk_bv_pair() : l_bytes(0), bit_bytes(nullptr){}
	~yon_blk_bv_pair(){ delete [] this->bit_bytes; }

	void clear(void){
		this->pattern.clear();
		this->l_bytes = 0;
		delete [] this->bit_bytes;
		this->bit_bytes = nullptr;
	}

	yon_blk_bv_pair& operator=(const yon_blk_bv_pair& other){
		delete [] this->bit_bytes;
		this->pattern   = other.pattern;
		this->l_bytes   = other.l_bytes;
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
		other.pattern.clear(); // Clear the src pattern vector.
		this->l_bytes = other.l_bytes;
		other.l_bytes = 0; // Clear the src byte length.
		return(*this);
	}

	/**<
	 * Predicate for a target local idx field in this pattern. Returns
	 * TRUE if the bit is set or FALSE otherwise. This function performs
	 * no checks for the target bit position being out-of-bounds.
	 * @param position Target bit position.
	 * @return         Returns TRUE if the bit is set or FALSE otherwise.
	 */
	inline bool operator[](const uint32_t position) const{ return((this->bit_bytes[position / 8] & (1 << (position % 8))) >> (position % 8)); }

	/**<
	 * Construct the lookup bit-vector for this object. This function needs
	 * to know the total number of fields that are set in the parent
	 * VariantBlockFooter structure as this will determine that byte-width
	 * of the bit-vector. Additionally, this function needs to be given a
	 * pointer to the map from global idx to local idx as the bit-vector
	 * bits corresponds to local idx predicates.
	 *
	 * Internally the bit-vector for all objects in a VariantBlockFooter
	 * structure has ceil(n_total_fields/8) bytes allocated for the base
	 * array.
	 * @param n_footer_total_fields Total number of fields set in the parent VariantBlockFooter.
	 * @param local_map             Pointer to map from global idx to local idx.
	 */
	void Build(const uint32_t n_footer_total_fields,
	           const std::unordered_map<uint32_t, uint32_t>* local_map)
	{
		if(this->pattern.size() == 0) return;
		assert(local_map != nullptr);

		// Determine the required byte width of the bit-vector.
		uint8_t bitvector_width = ceil((float)(n_footer_total_fields + 1) / 8);

		// Allocate new bit-vectors.
		delete [] this->bit_bytes;
		this->l_bytes = bitvector_width;
		this->bit_bytes = new uint8_t[bitvector_width];
		memset(this->bit_bytes, 0, sizeof(uint8_t)*bitvector_width);

		// Cycle over global idx values in the vector.
		for(uint32_t i = 0; i < this->pattern.size(); ++i){
			std::unordered_map<uint32_t, uint32_t>::const_iterator it = local_map->find(this->pattern[i]);
			assert(it != local_map->end());

			// Map from absolute key to local key.
			uint32_t local_key = it->second;
			assert(local_key <= n_footer_total_fields);

			// Set the target bit to TRUE at the local key position.
			this->bit_bytes[local_key/8] |= 1 << (local_key % 8);
		}
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const yon_blk_bv_pair& entry){
		io::SerializePrimitive(entry.l_bytes, buffer);
		buffer += (uint32_t)entry.pattern.size();
		for(uint32_t i = 0; i < entry.pattern.size(); ++i)
			io::SerializePrimitive(entry.pattern[i], buffer);

		for(uint32_t i = 0; i < entry.l_bytes; ++i)
			io::SerializePrimitive(entry.bit_bytes[i], buffer);


		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, yon_blk_bv_pair& entry){
		entry.pattern.clear();
		io::DeserializePrimitive(entry.l_bytes, buffer);
		uint32_t l_vector;
		buffer >> l_vector;
		//entry.pattern.resize(l_vector);
		for(uint32_t i = 0; i < l_vector; ++i){
			int temp;
			io::DeserializePrimitive(temp, buffer);
			//entry.pattern[i] = temp;
			entry.pattern.push_back(temp);
		}

		entry.bit_bytes = new uint8_t[entry.l_bytes];
		for(uint32_t i = 0; i < entry.l_bytes; ++i)
			io::DeserializePrimitive(entry.bit_bytes[i], buffer);

		return(buffer);
	}

public:
	std::vector<int> pattern; // vector of global idx values.
	uint8_t  l_bytes; // number of bytes in bit-vector.
	uint8_t* bit_bytes; // byte array interpreted as a bit-vector.
};

namespace containers {

struct VariantBlockFooter {
public:
	typedef VariantBlockFooter  self_type;
	typedef DataContainerHeader header_type;
	typedef std::unordered_map<uint32_t, uint32_t> map_type;
	typedef std::unordered_map<uint64_t, uint32_t> map_pattern_type;

public:
	VariantBlockFooter();
	~VariantBlockFooter();
	VariantBlockFooter(const self_type& other);
	VariantBlockFooter(self_type&& other) noexcept;
	VariantBlockFooter& operator=(const self_type& other);
	VariantBlockFooter& operator=(self_type&& other) noexcept;

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

}
}

#endif /* CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_ */
