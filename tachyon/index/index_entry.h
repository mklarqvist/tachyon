#ifndef INDEX_INDEX_ENTRY_H_
#define INDEX_INDEX_ENTRY_H_

#include <fstream>
#include <limits>

#include "../support/type_definitions.h"

namespace tachyon{
namespace index{

/**<
 * Sorted index entry
 */
struct IndexEntry{
	typedef IndexEntry self_type;

public:
	IndexEntry() :
		blockID(0),
		contigID(-1),
		n_variants(0),
		byte_offset(0),
		byte_offset_end(0),
		minPosition(0),
		maxPosition(0),
		minBin(std::numeric_limits<S32>::max()),
		maxBin(0)
	{}
	~IndexEntry(){}

	void reset(void){
		this->blockID         = 0;
		this->contigID        = -1;
		this->n_variants      = 0;
		this->byte_offset     = 0;
		this->byte_offset_end = 0;
		this->minPosition     = 0;
		this->maxPosition     = 0;
		this->minBin          = std::numeric_limits<S32>::max();
		this->maxBin          = 0;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.blockID),         sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),        sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),      sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minBin),          sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.maxBin),          sizeof(S32));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.blockID),         sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID),        sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),      sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minBin),          sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.maxBin),          sizeof(S32));

		return(stream);
	}

public:
	// Global index
	U64 blockID;        // this data is duplicated
	S32 contigID;       // this data is duplicated
	U32 n_variants;     // this data is duplicated
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;// tellg() position in stream for start of record in Tomahawk file
	U64 minPosition;	// smallest bp position
	U64 maxPosition;	// largest  bp position
	S32 minBin;
	S32 maxBin;
};

}
}

#endif
