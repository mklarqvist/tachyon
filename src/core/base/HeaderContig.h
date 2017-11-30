#ifndef VCF_VCFHEADERCONTIG_H_
#define VCF_VCFHEADERCONTIG_H_

namespace Tachyon{
namespace Core{

struct HeaderContig{
	typedef HeaderContig self_type;

public:
	HeaderContig() : l_name(0), bp_length(0), blocks(0){}
	~HeaderContig(){}

	inline void operator++(void){ ++this->blocks; }
	inline void operator--(void){ --this->blocks; }
	template <class T> inline void operator+=(const T value){ this->blocks += value; }
	template <class T> inline void operator-=(const T value){ this->blocks -= value; }

	friend std::ostream& operator<<(std::ostream& out, const self_type& contig){
		out << contig.l_name << '\t' << contig.name << '\t' << contig.bp_length << '\t' << contig.blocks;
		return(out);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.l_name), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.bp_length), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.blocks), sizeof(U32));
		stream.write(&entry.name[0], entry.l_name);
		return(stream);
	}

public:
	U16 l_name;
	U64 bp_length;
	// keep track of how many blocks we've seen for this contig
	// used during import
	U32 blocks;
	std::string name;
};

}
}



#endif /* VCF_VCFHEADERCONTIG_H_ */
