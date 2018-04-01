#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

#include <fstream>

#include "../../io/vcf/VCFHeaderLine.h"

namespace tachyon{
namespace core{

enum TACHYON_VARIANT_HEADER_FIELD_TYPE{
	YON_VCF_HEADER_INTEGER,
	YON_VCF_HEADER_FLOAT,
	YON_VCF_HEADER_FLAG,
	YON_VCF_HEADER_CHARACTER,
	YON_VCF_HEADER_STRING
};


/**<
 * FORMAT/FILTER/INFO field entry
 */
struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry() :
		IDX(0),
		primitive_type(0)
	{}

	HeaderMapEntry(const std::string& id, const S32& idx) :
		IDX(idx),
		primitive_type(0),
		ID(id)
	{}

	HeaderMapEntry(const std::string& id, const S32& idx, const S32& primitive_type) :
		IDX(idx),
		primitive_type(primitive_type),
		ID(id)
	{}

	HeaderMapEntry(const std::string& id) :
		IDX(0),
		primitive_type(0),
		ID(id)
	{}

	HeaderMapEntry(const self_type& other) :
		IDX(other.IDX),
		primitive_type(other.primitive_type),
		ID(other.ID)
	{}

	HeaderMapEntry(self_type&& other) :
		IDX(other.IDX),
		primitive_type(other.primitive_type),
		ID(other.ID)
	{}

	HeaderMapEntry& operator=(const self_type& other){
		this->IDX = other.IDX;
		this->primitive_type = other.primitive_type;
		this->ID  = other.ID;
		return(*this);
	}

	HeaderMapEntry& operator=(HeaderMapEntry&& other){
		this->IDX = other.IDX;
		this->primitive_type = other.primitive_type;
		this->ID  = other.ID;
		return(*this);
	}

	~HeaderMapEntry(){}

	inline const bool operator<(const self_type& other) const{
		return(this->IDX < other.IDX);
	}

	inline const bool operator>(const self_type& other) const{
		return(!this->operator<(other));
	}

	inline const TACHYON_VARIANT_HEADER_FIELD_TYPE getType(void) const{ return(TACHYON_VARIANT_HEADER_FIELD_TYPE(this->primitive_type)); }

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		const U32 l_ID = entry.ID.size();
		stream.write(reinterpret_cast<const char*>(&l_ID),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry),      sizeof(BYTE));
		stream.write(reinterpret_cast<const char*>(&entry.IDX),  sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.primitive_type),  sizeof(BYTE));
		stream.write(&entry.ID[0], entry.ID.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_ID = 0;
		stream.read(reinterpret_cast<char*>(&l_ID),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry),      sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&entry.IDX),  sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.primitive_type),  sizeof(BYTE));
		entry.ID.resize(l_ID);
		stream.read(&entry.ID[0], l_ID);

		return(stream);
	}

public:
	S32         IDX;
	BYTE        primitive_type;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
