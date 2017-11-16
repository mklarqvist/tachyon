#ifndef PERMUTATIONMANAGER_H_
#define PERMUTATIONMANAGER_H_

namespace Tomahawk{
namespace Core{

// Manages the PPA array and the Occ function
// of the PPA array
class PermutationManager{
	typedef PermutationManager self_type;

public:
	PermutationManager(){}

	// Lookup
	const BYTE& operator[](const BYTE& p) const{ return(*reinterpret_cast<const BYTE* const>(&this->PPA[p * sizeof(BYTE)])); }
	const U16&  operator[](const U16& p)  const{ return(*reinterpret_cast<const U16*  const>(&this->PPA[p * sizeof(U16)]));  }
	const U32&  operator[](const U32& p)  const{ return(*reinterpret_cast<const U32*  const>(&this->PPA[p * sizeof(U32)]));  }
	const U64&  operator[](const U64& p)  const{ return(*reinterpret_cast<const U64*  const>(&this->PPA[p * sizeof(U64)]));  }

public:
	// Not written to disk
	BYTE stide_size; // this equals to sizeof(T)
	U32 n_samples;   // redundancy but convenient

	// Data
	char* PPA;

	// Occ
};

}
}

#endif /* PERMUTATIONMANAGER_H_ */
