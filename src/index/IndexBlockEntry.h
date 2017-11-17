#ifndef INDEX_INDEXBLOCKENTRY_H_
#define INDEX_INDEXBLOCKENTRY_H_

#include <fstream>
#include <bitset>

#include "../io/BasicBuffer.h"
#include "../algorithm/OpenHashTable.h"
#include "../core/StreamContainer.h"
#include "../core/HashContainer.h"
#include "../core/base/EntryHotMeta.h"

namespace Tomahawk{
namespace Index{

// Controller
struct IndexBlockEntryController{
	typedef IndexBlockEntryController self_type;

public:
	IndexBlockEntryController():
		hasGT(0),
		isDiploid(0),
		unused(0)
	{}
	~IndexBlockEntryController(){}

	inline void clear(){ memset(this, 0, sizeof(BYTE)); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(U16));
		return(stream);
	}

public:
	U16 hasGT: 1,
	    isDiploid: 1,
		unused: 14;
};

// Size of entries in these records are
// inferred from the number of INFO/FORMAT/FILTER
// entries in all the records in a block
struct IndexBlockEntryBitvector{
	typedef IndexBlockEntryBitvector self_type;

public:
	IndexBlockEntryBitvector() : bit_bytes(nullptr){}
	~IndexBlockEntryBitvector(){ delete [] this->bit_bytes; }

	inline void update(const BYTE& value, const U32& pos){ this->bit_bytes[pos] = value; }
	inline void allocate(const U32& w){ this->bit_bytes = new BYTE[w]; memset(this->bit_bytes, 0, w); }

	template <class T>
	const bool operator[](const T& p) const{
		std::cerr << (U32)p << " byte: " << (U32)p/8 << " remainder: " << (U32)p%8 << " shift: " << std::bitset<8>((1 << (p % 8))) << std::endl;
		return((this->bit_bytes[p / 8] & (1 << (p % 8))) >> (p % 8));
	}

public:
	BYTE* bit_bytes;
};

struct IndexBlockEntryOffsets{
	typedef IndexBlockEntryOffsets self_type;
	typedef Core::StreamContainerHeader header_type;
	typedef Core::StreamContainerHeaderStride header_stride_type;

public:
	IndexBlockEntryOffsets(void) : key(0), header_stride(nullptr){}
	IndexBlockEntryOffsets(const U32& key, const header_type& h) : key(key), header(h), header_stride(nullptr){}
	IndexBlockEntryOffsets(const U32& key, const header_type& h, const header_stride_type& s) : key(key), header(h), header_stride(new header_stride_type){
		*this->header_stride = s;
	}
	~IndexBlockEntryOffsets(void){ delete this->header_stride; }

	bool update(const U32& key, const header_type& h){
		this->key = key;
		this->header = h;
		return true;
	}

	bool update(const U32& key, const header_type& h, const header_stride_type& s){
		this->key = key;
		this->header = h;
		delete this->header_stride;
		this->header_stride = new header_stride_type;
		*this->header_stride = s;
		return true;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.key), sizeof(U32));
		stream << entry.header;
		if(entry.header.controller.mixedStride)
			stream << *entry.header_stride;

		return(stream);
	}

public:
	U32 key;
	header_type header;
	header_stride_type* header_stride;
};

struct IndexBlockEntryBase{
	typedef IndexBlockEntryBase self_type;
	typedef IndexBlockEntryController controller_type;

public:
	IndexBlockEntryBase();
	virtual ~IndexBlockEntryBase();

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.contigID),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),  sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.offset_streams_begin),  sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_ppa), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_meta), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_meta_complex), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_rle), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_simple), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_info_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_streams), sizeof(U16));

		stream.write(reinterpret_cast<const char*>(&entry.n_info_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_patterns), sizeof(U16));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.contigID),  sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.offset_streams_begin), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_ppa), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_meta), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_meta_complex), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_rle), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_simple), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_info_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_streams), sizeof(U16));

		stream.read(reinterpret_cast<char*>(&entry.n_info_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_patterns), sizeof(U16));

		return(stream);
	}

	void reset(void){
		this->controller.clear();
		this->contigID = -1;
		this->n_variants = 0;
		this->offset_streams_begin = 0;
		this->l_ppa = 0;
		this->l_meta = 0;
		this->l_gt_rle = 0;
		this->l_gt_simple = 0;
		this->l_meta_complex = 0;
		this->n_info_streams = 0;
		this->n_format_streams = 0;
		this->n_filter_streams = 0;
		this->n_info_patterns = 0;
		this->n_format_patterns = 0;
		this->n_filter_patterns = 0;
	}

public:
	// 2 + 4 +2 + 4 + 5*4 + 3*2 = 38
	controller_type controller;    // bit flags
	S32 contigID;       // contig identifier
	U16 n_variants;     // number of variants in this block
	U32 offset_streams_begin; // start of PPA

	// Virtual offsets to the start of various
	// basic fields:
	// PPA, META, META_COMPLEX, GT_RLE, GT_SIMPLE
	U32 l_ppa;
	U32 l_meta;
	U32 l_meta_complex;
	U32 l_gt_rle;
	U32 l_gt_simple;

	// Number of INFO/FORMAT/FILTER streams
	// in this block
	//
	// Bit vectors + the stream keys = the vector of identifiers
	U16 n_info_streams;
	U16 n_format_streams;
	U16 n_filter_streams;

	// How many patterns is there?
	U16 n_info_patterns;
	U16 n_format_patterns;
	U16 n_filter_patterns;

	// END OF FIXED SIZE
};

struct IndexBlockEntry : public IndexBlockEntryBase{
private:
	typedef IndexBlockEntry self_type;
	typedef IndexBlockEntryBase base_type;
	typedef IndexBlockEntryController controller_type;
	typedef IndexBlockEntryBitvector bit_vector;
	typedef Hash::HashTable<U32, U32> hash_table;
	typedef std::vector<U32> id_vector;
	typedef std::vector< id_vector > pattern_vector;
	typedef Core::Support::HashContainer hash_container_type;
	typedef Core::Support::HashVectorContainer hash_vector_container_type;
	typedef IndexBlockEntryOffsets offset_type;

public:
	enum INDEX_BLOCK_TARGET{INDEX_INFO, INDEX_FORMAT, INDEX_FILTER};

public:
	IndexBlockEntry();
	~IndexBlockEntry();
	void reset(void);

	// Allocate offset vectors
	void allocateInfoOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->info_offsets;
		this->info_offsets = new offset_type[size];
	}

	void allocateFormatOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->format_offsets;
		this->format_offsets = new offset_type[size];
	}

	void allocateFilterOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->filter_offsets;
		this->filter_offsets = new offset_type[size];
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const IndexBlockEntryBase* const base = reinterpret_cast<const IndexBlockEntryBase* const>(&entry);
		stream << base;

		std::cerr << "Writing: " << entry.n_info_streams << " INFO headers..." << std::endl;
		for(U32 i = 0; i < entry.n_info_streams; ++i)
			stream << entry.info_offsets[i];

		std::cerr << "Writing: " << entry.n_format_streams << " FORMAT headers..." << std::endl;
		for(U32 i = 0; i < entry.n_format_streams; ++i)
			stream << entry.format_offsets[i];

		std::cerr << "Writing: " << entry.n_filter_streams << " FILTER headers..." << std::endl;
		for(U32 i = 0; i < entry.n_filter_streams; ++i)
			stream << entry.filter_offsets[i];

		// write
		if(entry.n_info_patterns > 0){
			const BYTE info_bitvector_width = ceil((float)entry.n_info_streams/8);
			std::cerr << "Writing info bit vectors... " << entry.n_info_patterns << std::endl;
			for(U32 i = 0; i < entry.n_info_patterns; ++i)
				stream.write((const char*)entry.info_bit_vectors[i].bit_bytes, info_bitvector_width);
		}

		if(entry.n_format_patterns > 0){
			std::cerr << "Writing format bit vectors... " << entry.n_format_patterns << std::endl;
			const BYTE format_bitvector_width = ceil((float)entry.n_format_streams/8);
			for(U32 i = 0; i < entry.n_format_patterns; ++i)
				stream.write((const char*)entry.format_bit_vectors[i].bit_bytes, format_bitvector_width);
		}


		if(entry.n_filter_patterns > 0){
			std::cerr << "Writing filter bit vectors... " << entry.n_filter_patterns << std::endl;
			const BYTE filter_bitvector_width = ceil((float)entry.n_filter_streams/8);
			for(U32 i = 0; i < entry.n_filter_patterns; ++i)
				stream.write((const char*)entry.filter_bit_vectors[i].bit_bytes, filter_bitvector_width);
		}


		return(stream);
	}

	/////////////////////////
	// Import functionality
	/////////////////////////
	// During import we need to intialize and
	// resize these these pointers to fit the
	// data we want to store
	bool constructBitVector(const INDEX_BLOCK_TARGET& target, hash_container_type& values, hash_vector_container_type& patterns);

private:
	bool __constructBitVector(bit_vector*& target, hash_container_type& values, hash_vector_container_type& patterns);

public:
	// Virtual offsets into various
	// INFO/FORMAT/FILTER streams
	//
	// This mean local key is implicit
	// GLOBAL KEY | OFFSET
	// These contain the local map to a stream ID
	// e.g. 15 -> 0, 18 -> 1, 8 -> 2 etc.
	offset_type* info_offsets;
	offset_type* format_offsets;
	offset_type* filter_offsets;

	// Structure of bit-vectors
	bit_vector* info_bit_vectors;
	bit_vector* format_bit_vectors;
	bit_vector* filter_bit_vectors;
};

}
}

#endif /* INDEX_IndexBLOCKENTRY_H_ */
