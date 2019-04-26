#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "genotypes.h"
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
	GenotypeSorter(const uint64_t n_samples);
	~GenotypeSorter();

	/**<
	 * Restore the initial state without releasing memory. This function
	 * does not need to be called every iteration as data is overwritten.
	 */
	void reset(void);

	// Accessor
	inline const uint64_t& GetNumberSamples(void) const { return(this->n_samples); }

	/**<
	 * Allocates memory to allocate for a set number of samples.
	 * Deletes existing data without consideration.
	 * @param n_samples Number of samples to allocate memory for.
	 */
	void SetSamples(const uint64_t n_samples);

	/**<
	 * Permutes the genotypes in the provided VcfContainer. This
	 * functions needs to have knowledge of the global Vcf header
	 * to locate the Format:GT field identifier. Genotypes will be
	 * permuted according to their genotypes. This function respects
	 * the base ploidy.
	 * @param vcf_container Src VcfContainer harbouring the Format:GT entries.
	 * @param vcf_header    Src global Vcf header.
	 * @return              Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(const vcf_container_type& vcf_container, io::VcfHeader& vcf_header);

	/**<
	 * Permutes the genotypes in the provided array of yon1_vnt_t records.
	 * Genotypes will be permuted according to their genotypes. This function respects
	 * the base ploidy and phasing.
	 * @param rec   Src array of yon1_vnt_t records.
	 * @param n_rec Number of records provided in rec.
	 * @return      Returns TRUE upon success or FALSE otherwise.
	 */
	bool Build(const yon1_vnt_t* rec, const uint32_t n_rec);
	void Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa);

public:
	uint64_t      n_samples; // total number of entries in file
	yon_gt_ppa    permutation_array;
	yon_radix_gt* gt_pattern;
	uint8_t       gt_remap[256];
};

}
}

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
