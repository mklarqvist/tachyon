#ifndef ALGORITHM_PERMUTATIONMANAGER_H_
#define ALGORITHM_PERMUTATIONMANAGER_H_

namespace Tachyon{
namespace Core{

// Manages the PPA array and the Occ function
// of the PPA array
class PermutationManager{
	typedef PermutationManager self_type;

public:
	PermutationManager() : stride_size(sizeof(U32)), n_samples(0), PPA(nullptr){}
	PermutationManager(const U32 n_samples) : stride_size(sizeof(U32)), n_samples(n_samples), PPA(new char[sizeof(U32)*n_samples]){}
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
	// Convenience function used during import only
	U32* get(void) const{ return(reinterpret_cast<U32*>(this->PPA)); }
	U32& operator[](const U32& p) const{ return(*reinterpret_cast<U32*>(&this->PPA[p * sizeof(U32)])); }
	//bool validate(void);

public:
	// Not written to disk
	BYTE stride_size; // this equals to sizeof(T)
	U32 n_samples;   // redundancy but convenient

	// Data
	char* PPA;
};

}
}

#endif /* ALGORITHM_PERMUTATIONMANAGER_H_ */
