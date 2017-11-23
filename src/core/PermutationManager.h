#ifndef ALGORITHM_PERMUTATIONMANAGER_H_
#define ALGORITHM_PERMUTATIONMANAGER_H_

#include <fstream>

namespace Tachyon{
namespace Core{

// Manages the PPA array and the Occ function
// of the PPA array
class PermutationManager{
	typedef PermutationManager self_type;

public:
	PermutationManager() : stride_size(sizeof(U32)), n_samples(0), c_length(0), crc(0), PPA(nullptr){}
	PermutationManager(const U32 n_samples) : stride_size(sizeof(U32)), n_samples(n_samples), c_length(0), crc(0), PPA(new char[sizeof(U32)*n_samples]){}
	~PermutationManager(){ delete [] this->PPA; }

	void setSamples(const U32 n_samples){
		this->n_samples = n_samples;
		delete [] this->PPA;

		this->PPA = new char[sizeof(U32)*n_samples];
		for(U32 i = 0; i < this->n_samples; ++i)
			(*this)[i] = i;
	}

	void reset(void){
		for(U32 i = 0; i < this->n_samples; ++i)
			(*this)[i] = i;
	}

	// Lookup
	// Convenience function used during import
	U32* get(void) const{ return(reinterpret_cast<U32*>(this->PPA)); }
	U32& operator[](const U32& p) const{ return(*reinterpret_cast<U32*>(&this->PPA[p * sizeof(U32)])); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& manager){
		stream << manager.c_length;
		stream << manager.crc;
		stream.write(manager.PPA, sizeof(U32)*manager.n_samples);
		return(stream);
	}

public:
	// Not written to disk
	BYTE stride_size; // this equals to sizeof(T)
	U32 n_samples;   // redundancy but convenient

	// Data
	U32 c_length;
	// Uncompressed length is implicit
	// as being SAMPLES * sizeof(U32)
	// U32 u_length;
	U32 crc;
	char* PPA;
};

}
}

#endif /* ALGORITHM_PERMUTATIONMANAGER_H_ */
