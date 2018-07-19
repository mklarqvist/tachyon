#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

#include <fstream>

#include "support/enums.h"
#include "io/basic_buffer.h"
#include "io/vcf/vcf_header_line.h"

namespace tachyon{
namespace core{


/**<
 * FORMAT/FILTER/INFO field entry
 */
struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry();
	HeaderMapEntry(const std::string& id, const S32& idx);
	HeaderMapEntry(const std::string& id, const S32& idx, const S32& primitive_type);
	HeaderMapEntry(const std::string& id);
	HeaderMapEntry(const self_type& other);
	HeaderMapEntry(self_type&& other);
	HeaderMapEntry& operator=(const self_type& other);
	HeaderMapEntry& operator=(HeaderMapEntry&& other);
	~HeaderMapEntry() = default;

	inline bool operator<(const self_type& other) const{ return(this->IDX < other.IDX); }
	inline bool operator>(const self_type& other) const{ return(!this->operator<(other));}
	inline TACHYON_VARIANT_HEADER_FIELD_TYPE getType(void) const{ return(TACHYON_VARIANT_HEADER_FIELD_TYPE(this->primitive_type)); }

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		const U32 l_ID = entry.ID.size();
		stream.write(reinterpret_cast<const char*>(&l_ID),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.IDX),  sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.primitive_type),  sizeof(BYTE));
		stream.write(&entry.ID[0], entry.ID.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_ID = 0;
		stream.read(reinterpret_cast<char*>(&l_ID),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.IDX),  sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.primitive_type),  sizeof(BYTE));
		entry.ID.resize(l_ID);
		stream.read(&entry.ID[0], l_ID);

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry){
		buffer += (U32)entry.ID.size();
		buffer += entry.IDX;
		buffer += entry.primitive_type;
		buffer.Add(&entry.ID[0], entry.ID.size());
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& entry){
		U32 l_ID = 0;
		buffer >> l_ID;
		buffer >> entry.IDX;
		buffer >> entry.primitive_type;
		entry.ID.resize(l_ID);
		buffer.read(&entry.ID[0], l_ID);
		return(buffer);
	}

public:
	S32         IDX;
	BYTE        primitive_type;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
