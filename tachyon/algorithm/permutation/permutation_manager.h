#ifndef ALGORITHM_PERMUTATIONMANAGER_H_
#define ALGORITHM_PERMUTATIONMANAGER_H_

#include <fstream>

#include "../../io/basic_buffer.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"

#include "../../containers/components/datacontainer_header.h"

namespace tachyon{
namespace algorithm{

// Manages the PPA array and the Occ function
// of the PPA array
class PermutationManager{
	typedef PermutationManager self_type;
	typedef io::BasicBuffer    buffer_type;
	typedef containers::DataContainerHeader header_type;

public:
	PermutationManager();
	PermutationManager(const U32 n_samples);
	~PermutationManager();

	void setSamples(const U32 n_samples);
	void reset(void);
	bool generateCRC(void);
	inline const U32 getObjectSize(void) const{ return(sizeof(U64) + this->header.data_header.cLength); }

	// Lookup
	// Convenience function used during import
	inline U32* get(void){ return(reinterpret_cast<U32*>(this->PPA.buffer)); }
	inline const U32* get(void) const{ return(reinterpret_cast<U32*>(this->PPA.buffer)); }
	inline U32& operator[](const U32& p){ return(*reinterpret_cast<U32*>(&this->PPA.buffer[p * sizeof(U32)])); }
	inline const U32& operator[](const U32& p) const{ return(*reinterpret_cast<U32*>(&this->PPA.buffer[p * sizeof(U32)])); }


private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& manager){
		stream.write(reinterpret_cast<const char*>(&manager.n_samples),sizeof(U64));
		//stream << manager.header;
		//stream.write(reinterpret_cast<const char*>(&manager.u_length), sizeof(U32));
		//stream.write(reinterpret_cast<const char*>(&manager.c_length), sizeof(U32));
		//stream.write(reinterpret_cast<const char*>(&manager.crc),      sizeof(U32));
		stream.write(manager.PPA.buffer, manager.header.data_header.cLength);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& manager){
		stream.read(reinterpret_cast<char*>(&manager.n_samples),sizeof(U64));
		//stream >> manager.header;
		//stream.read(reinterpret_cast<char*>(&manager.u_length), sizeof(U32));
		//stream.read(reinterpret_cast<char*>(&manager.c_length), sizeof(U32));
		//stream.read(reinterpret_cast<char*>(&manager.crc),      sizeof(U32));
		manager.PPA.resize(manager.header.data_header.uLength);
		stream.read(manager.PPA.buffer, manager.header.data_header.cLength);
		manager.PPA.n_chars = manager.header.data_header.cLength;
		return(stream);
	}

public:
	U64 n_samples; // redundancy but convenient
	header_type header;
	//U32 u_length;
	//U32 c_length;
	//U32 crc;
	buffer_type PPA;
};

}
}

#endif /* ALGORITHM_PERMUTATIONMANAGER_H_ */
