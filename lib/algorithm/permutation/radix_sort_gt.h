#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "io/bcf/bcf_reader.h"
#include "permutation_manager.h"
#include "core/genotype_summary.h"

#include "io/htslib_integration.h"

namespace tachyon {
namespace algorithm {

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * bi-allelic diploid.
 */
class RadixSortGT {
	typedef RadixSortGT        self_type;
	typedef bcf::BCFReader     bcf_reader_type;
	typedef bcf::BCFEntry      bcf_entry_type;
	typedef PermutationManager manager_type;

	typedef containers::VcfContainer vcf_container_type;

public:
	RadixSortGT();
	RadixSortGT(const U64 n_samples);
	~RadixSortGT();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	void reset(void);
	void SetSamples(const U64 n_samples);

	bool Build(const vcf_container_type& vcf_container);

	// Construct given a reader with a block
	// of BCF entries loaded in it
	bool Build(const bcf_reader_type& reader);
	bool Update(const bcf_entry_type& entry);

	inline const U64& GetSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }

	/**<
	* Static function that calculates the 64-bit hash value for the target
	* FORMAT/FILTER/INFO vector of id fields. The id fields must be of type
	* int (S32). Example of using this function:
	*
	* const U64 hash_value = VariantImporter::HashIdentifiers(id_vector);
	*
	* @param id_vector Input vector of FORMAT/FILTER/INFO identifiers.
	* @return          Returns a 64-bit hash value.
	*/
	template <class T>
	static U64 HashGenotypes(const T* genotype_data, const U32 l_data){
		XXH64_state_t* const state = XXH64_createState();
		if (state==NULL) abort();

		XXH_errorcode const resetResult = XXH64_reset(state, BCF_HASH_SEED);
		if (resetResult == XXH_ERROR) abort();

		for(U32 i = 0; i < l_data; ++i){
			XXH_errorcode const addResult = XXH64_update(state, (const void*)genotype_data, sizeof(T));
			if (addResult == XXH_ERROR) abort();
		}

		U64 hash = XXH64_digest(state);
		XXH64_freeState(state);

		return hash;
	}

public:
	U64           n_samples; // total number of entries in file
	U32           position;  // number of entries parsed
	U32           p_i[9];    // number of entries in bin i
	BYTE*         GT_array;  // packed genotype array
	U32**         bins;      // bin i
	manager_type* manager;   // permutation manager
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
