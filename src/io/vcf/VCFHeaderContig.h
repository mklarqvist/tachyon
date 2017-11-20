#ifndef VCF_VCFHEADERCONTIG_H_
#define VCF_VCFHEADERCONTIG_H_

namespace Tachyon{
namespace VCF{

struct VCFHeaderContig{
	typedef VCFHeaderContig self_type;

public:
	VCFHeaderContig() : length(0), blocks(0){}
	~VCFHeaderContig(){}

	inline void operator++(void){ ++this->blocks; }
	inline void operator--(void){ --this->blocks; }
	template <class T> inline void operator+=(const T value){ this->blocks += value; }
	template <class T> inline void operator-=(const T value){ this->blocks -= value; }

	friend std::ostream& operator<<(std::ostream& out, const self_type& contig){
		out << contig.name << '\t' << contig.length;
		return(out);
	}

public:
	std::string name;
	U32 length;
	// keep track of how many blocks we've seen for this contig
	// used during import
	U32 blocks;
};

}
}



#endif /* VCF_VCFHEADERCONTIG_H_ */
