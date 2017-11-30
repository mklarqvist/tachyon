#ifndef CORE_BASE_HEADERMAPENTRY_H_
#define CORE_BASE_HEADERMAPENTRY_H_

namespace Tachyon{
namespace Core{

struct HeaderMapEntry{
private:
	typedef HeaderMapEntry self_type;

public:
	HeaderMapEntry() : TYPE(){}
	HeaderMapEntry(const std::string& id, const BYTE& type) : TYPE(type), ID(id){}
	~HeaderMapEntry(){}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const U32 l_ID = entry.ID.size();
		stream.write(reinterpret_cast<const char*>(&l_ID),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.TYPE), sizeof(BYTE));
		stream.write(&entry.ID[0], entry.ID.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_ID = 0;
		stream.read(reinterpret_cast<char*>(&l_ID),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.TYPE), sizeof(BYTE));
		entry.ID.resize(l_ID);
		stream.read(&entry.ID[0], l_ID);
		return(stream);
	}

public:
	BYTE TYPE;
	std::string ID;
};

}
}

#endif /* CORE_BASE_HEADERMAPENTRY_H_ */
