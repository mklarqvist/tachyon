#ifndef INDEX_INDEX_INDEX_ENTRY_H_
#define INDEX_INDEX_INDEX_ENTRY_H_

#include "index/index_entry.h"

namespace tachyon{
namespace index{

/**<
 * Index of index entries
 */
struct IndexIndexEntry{
private:
	typedef IndexIndexEntry self_type;
	typedef IndexEntry entry_type;

public:
	IndexIndexEntry() :
		contigID(-1),
		n_blocks(0),
		n_variants(0),
		byte_offset_begin(0),
		byte_offest_end(0),
		minPosition(0),
		maxPosition(0)
	{}

	IndexIndexEntry(const self_type& other) :
		contigID(other.contigID),
		n_blocks(other.n_blocks),
		n_variants(other.n_variants),
		byte_offset_begin(other.byte_offset_begin),
		byte_offest_end(other.byte_offest_end),
		minPosition(other.minPosition),
		maxPosition(other.maxPosition)
	{

	}

	~IndexIndexEntry(){}

	void reset(void){
		this->contigID          = -1;
		this->n_blocks          = 0;
		this->n_variants        = 0;
		this->byte_offset_begin = 0;
		this->byte_offest_end   = 0;
		this->minPosition       = 0;
		this->maxPosition       = 0;
	}

	void operator()(const entry_type& entry){
		this->contigID          = entry.contigID;
		this->n_blocks          = 1;
		this->n_variants        = entry.n_variants;
		this->byte_offset_begin = entry.byte_offset;
		this->byte_offest_end   = entry.byte_offset_end;
		this->minPosition       = entry.minPosition;
		this->maxPosition       = entry.maxPosition;
	}

	inline bool operator==(const entry_type& entry) const{ return(this->contigID == entry.contigID); }
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	bool operator==(const self_type& other) const{
		if(this->contigID          != other.contigID) return false;
		if(this->n_blocks          != other.n_blocks) return false;
		if(this->n_variants        != other.n_variants) return false;
		if(this->byte_offset_begin != other.byte_offset_begin) return false;
		if(this->byte_offest_end   != other.byte_offest_end) return false;
		if(this->minPosition       != other.minPosition) return false;
		if(this->maxPosition       != other.maxPosition) return false;
		return true;
	}

	void operator+=(const entry_type& entry){
		++this->n_blocks;
		this->n_variants     += entry.n_variants;
		this->byte_offest_end = entry.byte_offset_end;
		this->maxPosition     = entry.maxPosition;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.contigID),         sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.n_blocks),         sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_begin),sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offest_end),  sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),      sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),      sizeof(U64));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.contigID),         sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.n_blocks),         sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_begin),sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offest_end),  sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),      sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),      sizeof(U64));

		return(stream);
	}

public:
	S32 contigID;           // contig ID
	U32 n_blocks;           // number of index entries for this contigID
	U64 n_variants;         // total number of variants
	U64 byte_offset_begin;  // start virtual offset
	U64 byte_offest_end;    // end virtual offset
	U64 minPosition;        // smallest bp position observed
	U64 maxPosition;        // largest  bp position observed
};

}
}


#endif /* INDEX_INDEX_INDEX_ENTRY_H_ */
