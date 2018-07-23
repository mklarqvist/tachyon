#ifndef CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_
#define CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_

#include "data_block_bitvector.h"
#include "data_container_header.h"
#include "containers/hash_container.h"
#include "io/basic_buffer.h"

namespace tachyon {
namespace containers {

// It is possible of getting mapping local indices to global IDX
// for either  FILTER/FORMAT/INFO fields by iterating over the
// relevant DataContainerHeader structures and tracking the incremental
// position they occur at (local IDX) and read the global IDX in the
// header structure.
struct VariantBlockFooter {
private:
	typedef VariantBlockFooter        self_type;
	typedef DataBlockBitvector        bit_vector;
	typedef std::vector<U32>          id_vector;
	typedef std::vector< id_vector >  pattern_vector;
	typedef containers::HashContainer hash_container_type;
	typedef containers::HashVectorContainer hash_vector_container_type;
	typedef DataContainerHeader       header_type;

public:
	// Internal use only
	enum INDEX_BLOCK_TARGET{INDEX_INFO, INDEX_FORMAT, INDEX_FILTER};

public:
	VariantBlockFooter();
	~VariantBlockFooter();
	void reset(void);

	// Allocate offset vectors
	inline void AllocateInfoHeaders(const U32 n_info_streams){
		delete [] this->info_offsets;
		if(n_info_streams == 0){
			this->info_offsets = nullptr;
			return;
		}
		this->info_offsets = new header_type[n_info_streams];
	}

	inline void AllocateFormatHeaders(const U32 n_format_streams){
		delete [] this->format_offsets;
		if(n_format_streams == 0){
			this->format_offsets = nullptr;
			return;
		}
		this->format_offsets = new header_type[n_format_streams];
	}

	inline void AllocateFilterHeaders(const U32 n_filter_streams){
		delete [] this->filter_offsets;
		if(n_filter_streams == 0){
			this->filter_offsets = nullptr;
			return;
		}
		this->filter_offsets = new header_type[n_filter_streams];
	}

	/**<
	 * Wrapper function for allocating new offset objects for INFO,
	 * FORMAT, and FILTER patterns and streams
	 * @param n_info_streams   Number of unique info streams
	 * @param n_format_streams Number of unique format streams
	 * @param n_filter_streams Number of unique filter streams
	 */
	inline void AllocateHeaders(const U32 n_info_streams,
		                        const U32 n_format_streams,
		                        const U32 n_filter_streams)
	{
		this->AllocateInfoHeaders(n_info_streams);
		this->AllocateFormatHeaders(n_format_streams);
		this->AllocateFilterHeaders(n_filter_streams);
	}

	/**<
	 * Wrapper function for constructing INFO/FORMAT/FILTER pattern
	 * set-membership bit-vectors
	 * @param target   Target group (INFO/FORMAT/FILTER)
	 * @param values   Hash container of values
	 * @param patterns Hash container of patterns
	 * @return         Returns TRUE upon success or FALSE otherwise
	 */
	bool constructBitVector(const INDEX_BLOCK_TARGET& target,
	                        hash_container_type& values,
	                        hash_vector_container_type& patterns);

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry);
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& entry);

private:
	/**<
	 *
	 * @param target
	 * @param offset
	 * @param values
	 * @param patterns
	 * @return
	 */
	bool __constructBitVector(bit_vector*& target,
	                          header_type* offset,
	                          hash_container_type& values,
	                          hash_vector_container_type& patterns);

public:
	// Utility members. l_*_bitvector stores the byte-length
	// of each of the bit-vectors described below.
	// Not written or read from disk. Used internally only
	// Todo: These should be embedded directly into the bit-vector
	//       structure.
	BYTE l_info_bitvector;
	BYTE l_format_bitvector;
	BYTE l_filter_bitvector;

	// Critical values used to track the number of data streams
	// that is available for each possible type. The n_*_streams
	// fields corresponds to the number of different byte streams
	// that are set in this block. The n_*_patterns corresponds
	// to the number of unique vectors of field identifiers that
	// occurred in the block.
	U16 n_info_streams;
	U16 n_format_streams;
	U16 n_filter_streams;
	U16 n_info_patterns;
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
	header_type  offset_ppa;
	header_type  offset_meta_contig;
	header_type  offset_meta_position;
	header_type  offset_meta_refalt;
	header_type  offset_meta_controllers;
	header_type  offset_meta_quality;
	header_type  offset_meta_names;
	header_type  offset_meta_alleles;
	header_type  offset_meta_info_id;
	header_type  offset_meta_format_id;
	header_type  offset_meta_filter_id;
	header_type  offset_gt_8b;
	header_type  offset_gt_16b;
	header_type  offset_gt_32b;
	header_type  offset_gt_64b;
	header_type  offset_gt_simple8;
	header_type  offset_gt_simple16;
	header_type  offset_gt_simple32;
	header_type  offset_gt_simple64;
	header_type  offset_gt_helper;
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
	bit_vector*  info_bit_vectors;
	bit_vector*  format_bit_vectors;
	bit_vector*  filter_bit_vectors;
};

}
}



#endif /* CONTAINERS_COMPONENTS_VARIANT_BLOCK_FOOTER_H_ */
