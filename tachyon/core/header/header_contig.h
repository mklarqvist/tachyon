#ifndef CORE_HEADERCONTIG_H_
#define CORE_HEADERCONTIG_H_

#include <fstream>
#include "../../io/basic_buffer.h"

namespace tachyon{
namespace core{

struct HeaderContig{
	typedef HeaderContig self_type;

public:
	HeaderContig() : contigID(0), bp_length(0), n_blocks(0){}
	~HeaderContig(){}

	inline void operator++(void){ ++this->n_blocks; }
	inline void operator--(void){ --this->n_blocks; }
	template <class T> inline void operator+=(const T value){ this->n_blocks += value; }
	template <class T> inline void operator-=(const T value){ this->n_blocks -= value; }

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		const U32 l_name = entry.name.size();
		stream.write(reinterpret_cast<const char*>(&l_name),          sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),  sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.bp_length), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_blocks),  sizeof(U32));
		stream.write(&entry.name[0], entry.name.size());
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		U32 l_name = 0;
		stream.read(reinterpret_cast<char*>(&l_name),          sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.contigID),  sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.bp_length), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_blocks),  sizeof(U32));
		entry.name.resize(l_name);
		stream.read(&entry.name[0], l_name);

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& contig){
		buffer += (U32)contig.name.size();
		buffer += contig.contigID;
		buffer += contig.bp_length;
		buffer += contig.n_blocks;
		buffer.Add(&contig.name[0], contig.name.size());
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& contig){
		U32 l_name;
		buffer >> l_name;
		buffer >> contig.contigID;
		buffer >> contig.bp_length;
		buffer >> contig.n_blocks;
		contig.name.resize(l_name);
		buffer.read(&contig.name[0], l_name);
		return(buffer);
	}

public:
	U32 contigID;
	U64 bp_length;
	U32 n_blocks;
	std::string name;
};

}
}

#endif /* CORE_HEADERCONTIG_H_ */
