#ifndef INDEX_INDEXBLOCKENTRY_H_
#define INDEX_INDEXBLOCKENTRY_H_

#include <fstream>
#include <bitset>

#include "../hash_container.h"
#include "../../io/basic_buffer.h"
#include "../components/datablock_bitvector.h"
#include "../components/datacontainer_header.h"
#include "../components/datacontainer_header_controller.h"

namespace tachyon{
namespace containers{

/** @brief Controller flags for an IndexBlockEntry
 * This structure is for internal use only and describes
 * various internal states as flags.
 */
struct DataBlockHeaderController{
	typedef DataBlockHeaderController self_type;

public:
	DataBlockHeaderController():
		hasGT(0),
		hasGTPermuted(0),
		anyEncrypted(0),
		unused(0)
	{}
	~DataBlockHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const U16 c = controller.hasGT |
                      controller.hasGTPermuted << 1 |
                      controller.anyEncrypted  << 2 |
                      controller.unused        << 3;

		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(U16));
		return(stream);
	}

public:
	U16 hasGT:         1,  // This block has GT FORMAT data
		hasGTPermuted: 1,  // have the GT fields been permuted
		anyEncrypted:  1,  // any data encrypted
		unused:        13; // reserved for future use
};

/** @brief Fixed-sized components of an IndexBlockEntry
 * For internal use only. This is a subcomponent of fixed
 * size describing:
 * 1) Contig, minimum and maximum genomic coordinates
 * 2) Offsets into the containers
 * 3) Number of containers and ID patterns
 * 4) Controller flags
 */
struct DataBlockHeader{
private:
	typedef DataBlockHeader              self_type;
	typedef DataBlockHeaderController    controller_type;
	typedef DataContainerHeader          header_type;

public:
	DataBlockHeader();
	~DataBlockHeader();

	inline const U32& size(void) const{ return(this->n_variants); }
	inline const S32& getContigID(void) const{ return(this->contigID); }
	inline const S64& getMinPosition(void) const{ return(this->minPosition); }
	inline const S64& getMaxPosition(void) const{ return(this->maxPosition); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.l_offset_footer),   sizeof(U32));
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U32));

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.l_offset_footer),   sizeof(U32));
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->l_offset_footer    = 0;
		this->controller.clear();
		this->contigID           = -1;
		this->minPosition        = 0;
		this->maxPosition        = 0;
		this->n_variants         = 0;
	}

public:
	// allows jumping to the next block when streaming
	// over the file and not using the index
	// EOF marker is at this position - sizeof(EOF marker)
	U32 l_offset_footer;
	controller_type controller;
	S32 contigID;       // contig identifier
	S64 minPosition;    // minimum coordinate in this block
	S64 maxPosition;    // maximum coordinate in this block
	U32 n_variants;     // number of variants in this block
};

struct DataBlockFooter{
private:
	typedef DataBlockFooter           self_type;
	typedef DataBlockBitvector        bit_vector;
	typedef hash::HashTable<U32, U32> hash_table;
	typedef std::vector<U32>          id_vector;
	typedef std::vector< id_vector >  pattern_vector;
	typedef containers::HashContainer hash_container_type;
	typedef containers::HashVectorContainer hash_vector_container_type;
	typedef DataContainerHeader       header_type;

public:
	// Internal use only
	enum INDEX_BLOCK_TARGET{INDEX_INFO, INDEX_FORMAT, INDEX_FILTER};

public:
	DataBlockFooter();
	~DataBlockFooter();
	void reset(void);

	// Allocate offset vectors
	inline void allocateInfoDiskOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->info_offsets;
		this->info_offsets = new header_type[size];
	}

	inline void allocateFormatDiskOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->format_offsets;
		this->format_offsets = new header_type[size];
	}

	inline void allocateFilterDiskOffsets(const U32& size){
		if(size == 0) return;
		delete [] this->filter_offsets;
		this->filter_offsets = new header_type[size];
	}

	inline void allocateDiskOffsets(const U32& info, const U32& format, const U32& filter){
		this->allocateInfoDiskOffsets(info);
		this->allocateFormatDiskOffsets(format);
		this->allocateFilterDiskOffsets(filter);
	}

	// During import we need to initialise and
	// resize these these pointers to fit the
	// data we want to store
	bool constructBitVector(const INDEX_BLOCK_TARGET& target, hash_container_type& values, hash_vector_container_type& patterns);

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& entry){
		buffer += entry.n_info_streams;
		buffer += entry.n_format_streams;
		buffer += entry.n_filter_streams;
		buffer += entry.n_info_patterns;
		buffer += entry.n_format_patterns;
		buffer += entry.n_filter_patterns;
		buffer += entry.offset_ppa;
		buffer += entry.offset_meta_contig;
		buffer += entry.offset_meta_position;
		buffer += entry.offset_meta_refalt;
		buffer += entry.offset_meta_controllers;
		buffer += entry.offset_meta_quality;
		buffer += entry.offset_meta_names;
		buffer += entry.offset_meta_alleles;
		buffer += entry.offset_gt_8b;
		buffer += entry.offset_gt_16b;
		buffer += entry.offset_gt_32b;
		buffer += entry.offset_gt_64b;
		buffer += entry.offset_gt_simple8;
		buffer += entry.offset_gt_simple16;
		buffer += entry.offset_gt_simple32;
		buffer += entry.offset_gt_simple64;
		buffer += entry.offset_gt_helper;
		buffer += entry.offset_meta_info_id;
		buffer += entry.offset_meta_format_id;
		buffer += entry.offset_meta_filter_id;

		for(U32 i = 0; i < entry.n_info_streams; ++i)
			buffer += entry.info_offsets[i];

		for(U32 i = 0; i < entry.n_format_streams; ++i)
			buffer += entry.format_offsets[i];

		for(U32 i = 0; i < entry.n_filter_streams; ++i)
			buffer += entry.filter_offsets[i];

		if(entry.n_info_patterns > 0){
			const BYTE info_bitvector_width = ceil((float)entry.n_info_streams/8);
			for(U32 i = 0; i < entry.n_info_patterns; ++i){
				buffer += entry.info_bit_vectors[i];
				buffer.Add((const char* const)&entry.info_bit_vectors[i].bit_bytes[0], info_bitvector_width);
			}
		}

		if(entry.n_format_patterns > 0){
			const BYTE format_bitvector_width = ceil((float)entry.n_format_streams/8);
			for(U32 i = 0; i < entry.n_format_patterns; ++i){
				buffer += entry.format_bit_vectors[i];
				buffer.Add((const char* const)&entry.format_bit_vectors[i].bit_bytes[0], format_bitvector_width);
			}
		}

		if(entry.n_filter_patterns > 0){
			const BYTE filter_bitvector_width = ceil((float)entry.n_filter_streams/8);
			for(U32 i = 0; i < entry.n_filter_patterns; ++i){
				buffer += entry.filter_bit_vectors[i];
				buffer.Add((const char* const)&entry.filter_bit_vectors[i].bit_bytes[0], filter_bitvector_width);
			}
		}

		return(buffer);
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.n_info_streams),    sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_streams),  sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_streams),  sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_info_patterns),   sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_patterns), sizeof(U16));

		stream << entry.offset_ppa;
		stream << entry.offset_meta_contig;
		stream << entry.offset_meta_position;
		stream << entry.offset_meta_refalt;
		stream << entry.offset_meta_controllers;
		stream << entry.offset_meta_quality;
		stream << entry.offset_meta_names;
		stream << entry.offset_meta_alleles;
		stream << entry.offset_gt_8b;
		stream << entry.offset_gt_16b;
		stream << entry.offset_gt_32b;
		stream << entry.offset_gt_64b;
		stream << entry.offset_gt_simple8;
		stream << entry.offset_gt_simple16;
		stream << entry.offset_gt_simple32;
		stream << entry.offset_gt_simple64;
		stream << entry.offset_gt_helper;
		stream << entry.offset_meta_info_id;
		stream << entry.offset_meta_format_id;
		stream << entry.offset_meta_filter_id;

		for(U32 i = 0; i < entry.n_info_streams; ++i)
			stream << entry.info_offsets[i];

		for(U32 i = 0; i < entry.n_format_streams; ++i)
			stream << entry.format_offsets[i];

		for(U32 i = 0; i < entry.n_filter_streams; ++i)
			stream << entry.filter_offsets[i];

		// write
		if(entry.n_info_patterns > 0){
			const BYTE info_bitvector_width = ceil((float)entry.n_info_streams/8);
			for(U32 i = 0; i < entry.n_info_patterns; ++i){
				stream << entry.info_bit_vectors[i];
				stream.write((const char*)entry.info_bit_vectors[i].bit_bytes, info_bitvector_width);
			}
		}

		if(entry.n_format_patterns > 0){
			const BYTE format_bitvector_width = ceil((float)entry.n_format_streams/8);
			for(U32 i = 0; i < entry.n_format_patterns; ++i){
				stream << entry.format_bit_vectors[i];
				stream.write((const char*)entry.format_bit_vectors[i].bit_bytes, format_bitvector_width);
			}
		}

		if(entry.n_filter_patterns > 0){
			const BYTE filter_bitvector_width = ceil((float)entry.n_filter_streams/8);
			for(U32 i = 0; i < entry.n_filter_patterns; ++i){
				stream << entry.filter_bit_vectors[i];
				stream.write((const char*)entry.filter_bit_vectors[i].bit_bytes, filter_bitvector_width);
			}
		}

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.n_info_streams),    sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_streams),  sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_streams),  sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_info_patterns),   sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_patterns), sizeof(U16));

		entry.l_info_bitvector   = ceil((float)entry.n_info_streams/8);
		entry.l_format_bitvector = ceil((float)entry.n_format_streams/8);
		entry.l_filter_bitvector = ceil((float)entry.n_filter_streams/8);

		stream >> entry.offset_ppa;
		stream >> entry.offset_meta_contig;
		stream >> entry.offset_meta_position;
		stream >> entry.offset_meta_refalt;
		stream >> entry.offset_meta_controllers;
		stream >> entry.offset_meta_quality;
		stream >> entry.offset_meta_names;
		stream >> entry.offset_meta_alleles;
		stream >> entry.offset_gt_8b;
		stream >> entry.offset_gt_16b;
		stream >> entry.offset_gt_32b;
		stream >> entry.offset_gt_64b;
		stream >> entry.offset_gt_simple8;
		stream >> entry.offset_gt_simple16;
		stream >> entry.offset_gt_simple32;
		stream >> entry.offset_gt_simple64;
		stream >> entry.offset_gt_helper;
		stream >> entry.offset_meta_info_id;
		stream >> entry.offset_meta_format_id;
		stream >> entry.offset_meta_filter_id;

		entry.info_offsets = new header_type[entry.n_info_streams];
		for(U32 i = 0; i < entry.n_info_streams; ++i)
			stream >> entry.info_offsets[i];

		entry.format_offsets = new header_type[entry.n_format_streams];
		for(U32 i = 0; i < entry.n_format_streams; ++i)
			stream >> entry.format_offsets[i];

		entry.filter_offsets = new header_type[entry.n_filter_streams];
		for(U32 i = 0; i < entry.n_filter_streams; ++i)
			stream >> entry.filter_offsets[i];

		if(entry.n_info_patterns > 0){
			const BYTE info_bitvector_width = ceil((float)entry.n_info_streams/8);
			entry.info_bit_vectors = new bit_vector[entry.n_info_patterns];
			for(U32 i = 0; i < entry.n_info_patterns; ++i){
				stream >> entry.info_bit_vectors[i];
				entry.info_bit_vectors[i].allocate(info_bitvector_width);
				stream.read((char*)entry.info_bit_vectors[i].bit_bytes, info_bitvector_width);

			}
		}

		if(entry.n_format_patterns > 0){
			const BYTE format_bitvector_width = ceil((float)entry.n_format_streams/8);
			entry.format_bit_vectors = new bit_vector[entry.n_format_patterns];
			for(U32 i = 0; i < entry.n_format_patterns; ++i){
				stream >> entry.format_bit_vectors[i];
				entry.format_bit_vectors[i].allocate(format_bitvector_width);
				stream.read((char*)entry.format_bit_vectors[i].bit_bytes, format_bitvector_width);
			}
		}

		if(entry.n_filter_patterns > 0){
			const BYTE filter_bitvector_width = ceil((float)entry.n_filter_streams/8);
			entry.filter_bit_vectors = new bit_vector[entry.n_filter_patterns];
			for(U32 i = 0; i < entry.n_filter_patterns; ++i){
				stream >> entry.filter_bit_vectors[i];
				entry.filter_bit_vectors[i].allocate(filter_bitvector_width);
				stream.read((char*)entry.filter_bit_vectors[i].bit_bytes, filter_bitvector_width);
			}
		}

		return(stream);
	}

private:
	/**<
	 *
	 * @param target
	 * @param offset
	 * @param values
	 * @param patterns
	 * @return
	 */
	bool __constructBitVector(bit_vector*& target, header_type* offset, hash_container_type& values, hash_vector_container_type& patterns);

public:
	U16 n_info_streams;
	U16 n_format_streams;
	U16 n_filter_streams;
	U16 n_info_patterns;
	U16 n_format_patterns;
	U16 n_filter_patterns;

	// Not written or read from disk
	// Used internally only
	BYTE l_info_bitvector;
	BYTE l_format_bitvector;
	BYTE l_filter_bitvector;

	// Headers of the various containers
	header_type  offset_ppa;
	header_type  offset_meta_contig;
	header_type  offset_meta_position;
	header_type  offset_meta_refalt;
	header_type  offset_meta_controllers;
	header_type  offset_meta_quality;
	header_type  offset_meta_names;
	header_type  offset_meta_alleles;
	header_type  offset_gt_8b;
	header_type  offset_gt_16b;
	header_type  offset_gt_32b;
	header_type  offset_gt_64b;
	header_type  offset_gt_simple8;
	header_type  offset_gt_simple16;
	header_type  offset_gt_simple32;
	header_type  offset_gt_simple64;
	header_type  offset_gt_helper;
	header_type  offset_meta_info_id;
	header_type  offset_meta_format_id;
	header_type  offset_meta_filter_id;
	header_type* info_offsets;
	header_type* format_offsets;
	header_type* filter_offsets;

	// Bit vectors
	bit_vector*  info_bit_vectors;
	bit_vector*  format_bit_vectors;
	bit_vector*  filter_bit_vectors;
};

}
}

#endif /* INDEX_IndexBLOCKENTRY_H_ */
