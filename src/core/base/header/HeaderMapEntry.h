#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

#include <fstream>

namespace Tachyon{
namespace Core{

struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry() : TYPE(0), IDX(0){}
	HeaderMapEntry(const std::string& id, const BYTE& type) : TYPE(type), IDX(0), ID(id){}
	HeaderMapEntry(const std::string& id, const BYTE& type, const S32& idx) : TYPE(type), IDX(idx), ID(id){}
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
		stream.write(reinterpret_cast<const char*>(&entry.TYPE), sizeof(BYTE));
		stream.write(reinterpret_cast<const char*>(&entry.IDX), sizeof(S32));
		stream.write(&entry.ID[0], entry.ID.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_ID = 0;
		stream.read(reinterpret_cast<char*>(&l_ID),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.TYPE), sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&entry.IDX), sizeof(S32));
		entry.ID.resize(l_ID);
		stream.read(&entry.ID[0], l_ID);

		return(stream);
	}

public:
	BYTE TYPE;
	S32 IDX;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
