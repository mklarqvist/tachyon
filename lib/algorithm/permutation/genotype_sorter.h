#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "core/genotypes.h"
#include "containers/vcf_container.h"

namespace tachyon {
namespace algorithm {

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * bi-allelic diploid.
 */
class GenotypeSorter {
public:
	typedef GenotypeSorter           self_type;
	typedef containers::VcfContainer vcf_container_type;

public:
	GenotypeSorter();
	GenotypeSorter(const U64 n_samples);
	~GenotypeSorter();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	void reset(void);
	void SetSamples(const U64 n_samples);

	bool Build(const vcf_container_type& vcf_container, io::VcfHeader& vcf_header);
	void Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa);

	inline const U64& GetNumberSamples(void) const{ return(this->n_samples); }

public:
	U64           n_samples; // total number of entries in file
	yon_gt_ppa    permutation_array;
	yon_radix_gt* gt_pattern;
	uint8_t gt_remap[256];
};

}
}

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
