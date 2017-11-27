#ifndef ALGORITHM_PERMUTATIONMANAGER_H_
#define ALGORITHM_PERMUTATIONMANAGER_H_

#include <fstream>

#include "../io/BasicBuffer.h"
#include "../third_party/zlib/zconf.h"
#include "../third_party/zlib/zlib.h"

namespace Tachyon{
namespace Core{

// Manages the PPA array and the Occ function
// of the PPA array
class PermutationManager{
	typedef PermutationManager self_type;
	typedef IO::BasicBuffer buffer_type;

public:
	PermutationManager() : stride_size(sizeof(U32)), n_samples(0), u_length(0), c_length(0), crc(0){}
	PermutationManager(const U32 n_samples) : stride_size(sizeof(U32)), n_samples(n_samples), u_length(0), c_length(0), crc(0), PPA(sizeof(U32)*n_samples){}
	~PermutationManager(){ this->PPA.deleteAll(); }

	void setSamples(const U32 n_samples){
		this->n_samples = n_samples;
		this->PPA.reset();
		this->PPA.resize(sizeof(S32)*n_samples);

		for(U32 i = 0; i < this->n_samples; ++i)
			this->PPA += (U32)i;
	}

	void reset(void){
		for(U32 i = 0; i < this->n_samples; ++i)
			(*this)[i] = i;

		this->u_length = 0;
		this->c_length = 0;
		this->crc = 0;
	}

	// Lookup
	// Convenience function used during import
	U32* get(void) const{ return(reinterpret_cast<U32*>(this->PPA.data)); }
	U32& operator[](const U32& p) const{ return(*reinterpret_cast<U32*>(&this->PPA.data[p * this->stride_size])); }

	bool generateCRC(void){
		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->PPA.data, this->n_samples*sizeof(U32));
		this->crc = crc;
		return true;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& manager){
		stream << manager.u_length;
		stream << manager.c_length;
		stream << manager.crc;
		stream.write(manager.PPA.data, manager.c_length);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& manager){
		stream >> manager.u_length;
		stream >> manager.c_length;
		stream >> manager.crc;
		manager.PPA.resize(manager.u_length);
		stream.read(manager.PPA.data, manager.c_length);
		return(stream);
	}

public:
	// Not written to disk
	BYTE stride_size; // this equals to sizeof(T)
	U32 n_samples;   // redundancy but convenient

	// Data
	U32 u_length;
	U32 c_length;
	U32 crc;
	buffer_type PPA;
};

}
}

#endif /* ALGORITHM_PERMUTATIONMANAGER_H_ */
