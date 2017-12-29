#ifndef INDEX_INDEXENTRY_H_
#define INDEX_INDEXENTRY_H_

#include <fstream>

#include "../support/TypeDefinitions.h"

namespace Tachyon{
namespace Index{

/**<
 * Sorted index entry
 */
struct IndexEntry{
	typedef IndexEntry self_type;

public:
	IndexEntry() :
		contigID(-1),
		n_variants(0),
		byte_offset(0),
		byte_offset_end(0),
		minPosition(0),
		maxPosition(0)
	{}
	~IndexEntry(){}

	void reset(void){
		this->contigID        = -1;
		this->n_variants      = 0;
		this->byte_offset     = 0;
		this->byte_offset_end = 0;
		this->minPosition     = 0;
		this->maxPosition     = 0;
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.contigID),        sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),      sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),     sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),     sizeof(U64));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.contigID),        sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),      sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),     sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),     sizeof(U64));

		return(stream);
	}

public:
	// Global index
	S32 contigID;       // this data is duplicated
	U16 n_variants;     // this data is duplicated
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;// tellg() position in stream for start of record in Tomahawk file
	U64 minPosition;	// smallest bp position
	U64 maxPosition;	// largest  bp position
};

}
}

#endif
