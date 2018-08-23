#ifndef INDEX_INDEX_INDEX_ENTRY_H_
#define INDEX_INDEX_INDEX_ENTRY_H_

#include "index/index_entry.h"

namespace tachyon{
namespace index{

struct IndexIndexEntry{
public:
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
		stream.write(reinterpret_cast<const char*>(&entry.contigID),         sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_blocks),         sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),       sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_begin),sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offest_end),  sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),      sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),      sizeof(uint64_t));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.contigID),         sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_blocks),         sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),       sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_begin),sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offest_end),  sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),      sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),      sizeof(uint64_t));

		return(stream);
	}

public:
	int32_t contigID;            // contig ID
	uint32_t n_blocks;           // number of index entries for this contig
	uint64_t n_variants;         // total number of variants
	uint64_t byte_offset_begin;  // start virtual offset
	uint64_t byte_offest_end;    // end virtual offset
	uint64_t minPosition;        // smallest bp position observed
	uint64_t maxPosition;        // largest  bp position observed
};

}
}


#endif /* INDEX_INDEX_INDEX_ENTRY_H_ */
