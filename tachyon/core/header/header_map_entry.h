#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

#include <fstream>

#include "../../io/vcf/VCFHeaderLine.h"

namespace tachyon{
namespace core{

enum TACHYON_VCF_HEADER_LINE_TYPE{
	YON_VCF_HEADER_UNKNOWN,
	YON_VCF_HEADER_FORMAT,
	YON_VCF_HEADER_FILTER,
	YON_VCF_HEADER_INFO,
	YON_VCF_HEADER_CONTIG
};

/**<
 * FORMAT/FILTER/INFO field entry
 */
struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry() :
		isFormat(0),
		isFilter(0),
		isInfo(0),
		unused(0),
		IDX(0)
	{}

	HeaderMapEntry(const std::string& id, const S32& idx) :
		isFormat(0),
		isFilter(0),
		isInfo(0),
		unused(0),
		IDX(idx),
		ID(id)
	{}

	HeaderMapEntry(const std::string& id, const BYTE isFormat, const BYTE isFilter, const BYTE isInfo) :
		isFormat(isFormat),
		isFilter(isFilter),
		isInfo(isInfo),
		unused(0),
		IDX(0),
		ID(id)
	{}

	HeaderMapEntry(const std::string& id, const S32& idx, const BYTE isFormat, const BYTE isFilter, const BYTE isInfo) :
		isFormat(isFormat),
		isFilter(isFilter),
		isInfo(isInfo),
		unused(0),
		IDX(idx),
		ID(id)
	{}

	HeaderMapEntry(const self_type& other) :
		isFormat(other.isFormat),
		isFilter(other.isFilter),
		isInfo(other.isInfo),
		unused(other.unused),
		IDX(other.IDX),
		ID(other.ID)
	{}

	HeaderMapEntry(self_type&& other) :
		isFormat(other.isFormat),
		isFilter(other.isFilter),
		isInfo(other.isInfo),
		unused(other.unused),
		IDX(other.IDX),
		ID(other.ID)
	{}

	HeaderMapEntry& operator=(const self_type& other){
		this->isFormat = other.isFormat;
		this->isFilter = other.isFilter;
		this->isInfo = other.isInfo;
		this->unused = other.unused;
		this->IDX = other.IDX;
		this->ID = other.ID;
		return(*this);
	}

	HeaderMapEntry& operator=(HeaderMapEntry&& other){
		this->isFormat = other.isFormat;
		this->isFilter = other.isFilter;
		this->isInfo = other.isInfo;
		this->unused = other.unused;
		this->IDX = other.IDX;
		this->ID = other.ID;
		return(*this);
	}

	~HeaderMapEntry(){}

	inline const bool operator<(const self_type& other) const{
		return(this->IDX < other.IDX);
	}

	inline const bool operator>(const self_type& other) const{
		return(!this->operator <(other));
	}

	inline bool isSet(const tachyon::core::TACHYON_VCF_HEADER_LINE_TYPE& type) const{
		switch(type){
		case(tachyon::core::TACHYON_VCF_HEADER_LINE_TYPE::YON_VCF_HEADER_FILTER): return(this->isFilter);
		case(tachyon::core::TACHYON_VCF_HEADER_LINE_TYPE::YON_VCF_HEADER_INFO): return(this->isInfo);
		case(tachyon::core::TACHYON_VCF_HEADER_LINE_TYPE::YON_VCF_HEADER_FORMAT): return(this->isFormat);
		default: return false;
		}
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const U32 l_ID = entry.ID.size();
		stream.write(reinterpret_cast<const char*>(&l_ID),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry),      sizeof(BYTE));
		stream.write(reinterpret_cast<const char*>(&entry.IDX),  sizeof(S32));
		stream.write(&entry.ID[0], entry.ID.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_ID = 0;
		stream.read(reinterpret_cast<char*>(&l_ID),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry),      sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&entry.IDX),  sizeof(S32));
		entry.ID.resize(l_ID);
		stream.read(&entry.ID[0], l_ID);

		return(stream);
	}

public:
	BYTE isFormat: 1,
         isFilter: 1,
         isInfo:   1,
         unused:   5;
	S32         IDX;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
