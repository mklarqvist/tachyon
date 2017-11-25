#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "../../io/bcf/BCFReader.h"
#include "../../core/PermutationManager.h"

// Todo: test
// if #common variants over some threshold then sort on those only
// beause they will dominate the RLE cost

namespace Tachyon {
namespace Algorithm {

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * bi-allelic diploid.
 */
class RadixSortGT {
	typedef RadixSortGT self_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry  bcf_entry_type;
	typedef Core::PermutationManager manager_type;

public:
	RadixSortGT();
	RadixSortGT(const U64 n_samples);
	~RadixSortGT();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	void reset(void);
	void setSamples(const U64 n_samples);

	// Construct given a reader with a block
	// of BCF entries loaded in it
	bool build(const bcf_reader_type& reader);
	bool update(const bcf_entry_type& entry);

	inline const U64& getSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }
	//inline U32& operator[](const U32& p) const{return((*this->manager)[p]); }
	//inline const U32* getPPA(void) const{return(this->ppa); }

	// Debug functions
	void outputGT(const bcf_reader_type& reader);

public:
	U64 n_samples; // total number of entries in file
	U32 position;  // number of entries parsed
	U32 p_i[9];    // number of entries in bin i
	//U32* ppa;      // position prefix array
	BYTE* GT_array;// packed genotype array
	U32** bins;    // bin i
	manager_type* manager; // permutation manager
	U64 cumulative_AAC;
	U64 cumulative_total;
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
