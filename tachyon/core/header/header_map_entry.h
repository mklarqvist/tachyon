#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

#include <fstream>

#include "../../io/vcf/VCFHeaderLine.h"

namespace tachyon{
namespace core{


/**<
 * FORMAT/FILTER/INFO field entry
 */
struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry() :
		IDX(0)
	{}

	HeaderMapEntry(const std::string& id, const S32& idx) :
		IDX(idx),
		ID(id)
	{}

	HeaderMapEntry(const std::string& id) :
		IDX(0),
		ID(id)
	{}

	HeaderMapEntry(const self_type& other) :
		IDX(other.IDX),
		ID(other.ID)
	{}

	HeaderMapEntry(self_type&& other) :
		IDX(other.IDX),
		ID(other.ID)
	{}

	HeaderMapEntry& operator=(const self_type& other){
		this->IDX = other.IDX;
		this->ID  = other.ID;
		return(*this);
	}

	HeaderMapEntry& operator=(HeaderMapEntry&& other){
		this->IDX = other.IDX;
		this->ID  = other.ID;
		return(*this);
	}

	~HeaderMapEntry(){}

	inline const bool operator<(const self_type& other) const{
		return(this->IDX < other.IDX);
	}

	inline const bool operator>(const self_type& other) const{
		return(!this->operator <(other));
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
	S32         IDX;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
