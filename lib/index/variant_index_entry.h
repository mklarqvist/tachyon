#ifndef INDEX_INDEX_ENTRY_H_
#define INDEX_INDEX_ENTRY_H_

#include <fstream>
#include <limits>

namespace tachyon{
namespace index{

struct VariantIndexEntry{
public:
	typedef VariantIndexEntry self_type;

public:
	VariantIndexEntry() :
		blockID(0),
		contigID(-1),
		n_variants(0),
		byte_offset(0),
		byte_offset_end(0),
		minPosition(0),
		maxPosition(0),
		minBin(std::numeric_limits<int32_t>::max()),
		maxBin(0)
	{}

	bool operator!=(const self_type& other) const{
		if(this->blockID         != other.blockID)         return true;
		if(this->contigID        != other.contigID)        return true;
		if(this->n_variants      != other.n_variants)      return true;
		if(this->byte_offset     != other.byte_offset)     return true;
		if(this->byte_offset_end != other.byte_offset_end) return true;
		if(this->minPosition     != other.minPosition)     return true;
		if(this->maxPosition     != other.maxPosition)     return true;
		if(this->minBin          != other.minBin)          return true;
		if(this->maxBin          != other.maxBin)          return true;
		return false;
	}

	inline bool operator==(const self_type& other) const{ return(!(*this != other)); }

	bool operator<(const self_type& other) const{
		if(this->blockID  < other.blockID)  return true;
		if(this->blockID  > other.blockID)  return false;
		if(this->contigID < other.contigID) return true;
		if(this->contigID > other.contigID) return false;
		return true;
	}

	bool operator<=(const self_type& other) const{
		if(this->blockID  <= other.blockID)  return true;
		if(this->blockID  >  other.blockID)  return false;
		if(this->contigID <= other.contigID) return true;
		if(this->contigID >  other.contigID) return false;
		return true;
	}

	~VariantIndexEntry(){}

	void reset(void){
		this->blockID         = 0;
		this->contigID        = -1;
		this->n_variants      = 0;
		this->byte_offset     = 0;
		this->byte_offset_end = 0;
		this->minPosition     = 0;
		this->maxPosition     = 0;
		this->minBin          = std::numeric_limits<int32_t>::max();
		this->maxBin          = 0;
	}

	std::ostream& print(std::ostream& stream) const{
		stream << blockID << '\t' << contigID << '\t' << n_variants << "\tOffset: " << byte_offset << '-' << byte_offset_end << " Position: " << minPosition << '-' << maxPosition << " Bins: " << minBin << '-' << maxBin;
		return(stream);
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.blockID),         sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),        sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),      sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.minBin),          sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.maxBin),          sizeof(int32_t));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.blockID),         sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.contigID),        sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),      sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.minBin),          sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.maxBin),          sizeof(int32_t));

		return(stream);
	}

public:
	uint64_t blockID; // this data is duplicated
	int32_t contigID; // this data is duplicated
	uint32_t n_variants; // this data is duplicated
	uint64_t byte_offset; // tellg() position in stream for start of record in Tachyon file
	uint64_t byte_offset_end;// tellg() position in stream for start of record in Tachyon file
	uint64_t minPosition;   // smallest bp position
	uint64_t maxPosition;  // largest  bp position
	int32_t minBin;
	int32_t maxBin;
};

}
}

#endif
