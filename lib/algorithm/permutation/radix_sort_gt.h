#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "permutation_manager.h"
#include "core/genotype_summary.h"
#include "core/genotypes.h"
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
	void Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa);

	inline const U64& GetNumberSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }

public:
	U64           n_samples; // total number of entries in file
	U32           position;  // number of entries parsed
	U32           p_i[9];    // number of entries in bin i
	BYTE*         GT_array;  // packed genotype array
	U32**         bins;      // bin i
	manager_type* manager;   // permutation manager

	//
	yon_gt_ppa    permutation_array;
	yon_radix_gt* gt_pattern;
	uint8_t gt_remap[256];
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
