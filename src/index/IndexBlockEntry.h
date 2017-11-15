#ifndef INDEX_INDEXBLOCKENTRY_H_
#define INDEX_INDEXBLOCKENTRY_H_

#include "../support/TypeDefinitions.h"

namespace Tomahawk{
namespace Index{

// Controller
struct IndexBlockEntryController{
	typedef IndexBlockEntryController self_type;

public:
	IndexBlockEntryController():
		mixedStride(0),
		signedness(0),
		type(0),
		method(0)
	{}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const BYTE* c = reinterpret_cast<const BYTE* const>(&controller);
		stream.write(reinterpret_cast<const char*>(&c), sizeof(BYTE));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		BYTE* c = reinterpret_cast<BYTE*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(BYTE));
		return(stream);
	}

public:
	BYTE mixedStride: 1,
	     signedness: 1,
		 type: 3,
		 method: 3;
};

// Block offset (key,key)-pair
struct IndexBlockEntryOffset{
	typedef IndexBlockEntryOffset self_type;

public:
	IndexBlockEntryOffset() : global_key(0), offset(0){}
	IndexBlockEntryOffset(const U32& global_key, const U32& offset) : global_key(global_key), offset(offset){}
	~IndexBlockEntryOffset(){}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& record){
		stream.write(reinterpret_cast<const char*>(&record.global_key), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&record.offset),     sizeof(U32));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& record){
		stream.read(reinterpret_cast<char*>(&record.global_key),sizeof(U32));
		stream.read(reinterpret_cast<char*>(&record.offset),    sizeof(U32));

		return(stream);
	}

public:
	U32 global_key;
	U32 offset;
};

// Size of entries in these records are
// inferred from the number of INFO/FORMAT/FILTER
// entries in all the records in a block
struct IndexBlockEntryBitvector{
	typedef IndexBlockEntryBitvector self_type;

public:
	// bool construct(const std::vector<U32>& values, const U32& width);
	inline void update(const BYTE& value, const U32& pos){ this->bit_bytes[pos] = value; }

	template <class T>
	FORCE_INLINE const bool operator[](const T& p) const{ return((this->bit_bytes[p / 8] & (1 << (p % 8))) >> (p % 8)); }

public:
	BYTE* bit_bytes;
};

struct IndexBlockEntryBase{
	typedef IndexBlockEntryBase self_type;
	typedef IndexBlockEntryController controller_type;

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.contigID),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),  sizeof(U16));

		stream.write(reinterpret_cast<const char*>(&entry.l_meta), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_meta_complex), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_rle), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_simple), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_info_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_streams), sizeof(U16));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.contigID),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),       sizeof(U16));

		stream.read(reinterpret_cast<char*>(&entry.l_meta), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_meta_complex), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_rle), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_simple), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_info_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_streams), sizeof(U16));

		return(stream);
	}

	void reset(void){
		this->contigID = 0;
		this->n_variants = 0;
		this->l_meta = 0;
		this->l_gt_rle = 0;
		this->l_gt_simple = 0;
		this->l_meta_complex = 0;
		this->n_info_streams = 0;
		this->n_format_streams = 0;
		this->n_filter_streams = 0;
	}

public:
	controller_type controller;    // bit flags
	S32 contigID;       // contig identifier
	U16 n_variants;     // number of variants in this block
	U32 offset_streams_begin; // start of PPA

	// Virtual offsets to the start of various
	// basic fields:
	// PPA, META, META_COMPLEX, GT_RLE, GT_SIMPLE
	/*
	 Streams
	 +---+---+---+---+---+---+---+---+---+---+
	 | 1 | 4             |                 ~~~
	 +---+---+---+---+---+---+---+---+---+---+
	 ^   ^               ^
	 |   |               |
	 CNT U_SIZE          DATA
	*/
	U32 l_ppa;
	U32 l_meta;
	U32 l_meta_complex;
	U32 l_gt_rle;
	U32 l_gt_simple;

	// Number of INFO/FORMAT/FILTER streams
	// in this block
	U16 n_info_streams;
	U16 n_format_streams;
	U16 n_filter_streams;
	// END OF FIXED SIZE

};

struct IndexBlockEntry{
	typedef IndexBlockEntry self_type;
	typedef IndexBlockEntryBase base_type;
	typedef IndexBlockEntryController controller_type;
	typedef IndexBlockEntryOffset offset_type;
	typedef IndexBlockEntryBitvector bit_vector;

public:
	IndexBlockEntry(){}
	~IndexBlockEntry(){}

	void reset(void){
		this->base.reset();

		delete [] this->info_offsets;
		delete [] this->format_offsets;
		delete [] this->filter_offsets;
		this->info_offsets = nullptr;
		this->format_offsets = nullptr;
		this->filter_offsets = nullptr;

		delete [] this->info_bit_vectors;
		delete [] this->format_bit_vectors;
		delete [] this->filter_bit_vectors;
		this->info_bit_vectors = nullptr;
		this->format_bit_vectors = nullptr;
		this->filter_bit_vectors = nullptr;
	}

public:
	base_type base;

	// Virtual offsets into various
	// INFO/FORMAT/FILTER streams

	// This mean local key is implicit
	// GLOBAL KEY | OFFSET
	offset_type* info_offsets;
	offset_type* format_offsets;
	offset_type* filter_offsets;
	/*
	 Data starts at offset
	 +---+---+---+---+---+---+---+---+
	 | 1 | 1 | 4 | 4 | 4 | 1 | 4 | 4 |
	 +---+---+---+---+---+---+---+---+
	 ^   ^   ^   ^   ^   ^   ^   ^
	 |   |   |   |   |   |   |   |
	 CNT STRIDE  CLENGTH CNT CLENGTH
	         OFFSET  ULENGTH     ULENGTH
	*/

	/*
	//Base struct
	BYTE controller;
	BYTE stride;
	U32 offset;
	U32 cLength;
	U32 uLength;
	// If extended
	BYTE controller2;
	U32 cLengthStride;
	U32 uLengthStride;
	*/

	// Structure of bit-vectors
	bit_vector* info_bit_vectors;
	bit_vector* format_bit_vectors;
	bit_vector* filter_bit_vectors;

	// Store associative maps
	// and bit vectors
	//
	// Remainder is
	// INFO and FORMAT and FILTER
	// fields
};

}
}

#endif /* INDEX_IndexBLOCKENTRY_H_ */
