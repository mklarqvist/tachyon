#ifndef CORE_HEADERCONTIG_H_
#define CORE_HEADERCONTIG_H_

namespace Tachyon{
namespace Core{

struct HeaderContig{
	typedef HeaderContig self_type;

public:
	HeaderContig() : bp_length(0), blocks(0){}
	~HeaderContig(){}

	inline void operator++(void){ ++this->blocks; }
	inline void operator--(void){ --this->blocks; }
	template <class T> inline void operator+=(const T value){ this->blocks += value; }
	template <class T> inline void operator-=(const T value){ this->blocks -= value; }

	friend std::ostream& operator<<(std::ostream& out, const self_type& contig){
		out << contig.name << '\t' << contig.bp_length << '\t' << contig.blocks;
		return(out);
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const U32 l_name = entry.name.size();
		stream.write(reinterpret_cast<const char*>(&l_name), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.bp_length), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.blocks), sizeof(U32));
		stream.write(&entry.name[0], entry.name.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_name = 0;
		stream.read(reinterpret_cast<char*>(&l_name),          sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.bp_length), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.blocks),    sizeof(U32));
		entry.name.resize(l_name);
		stream.read(&entry.name[0], l_name);

		return(stream);
	}

public:
	U64 bp_length;
	// keep track of how many blocks we've seen for this contig
	// used during import
	U32 blocks;
	std::string name;
};

}
}

#endif /* CORE_HEADERCONTIG_H_ */
