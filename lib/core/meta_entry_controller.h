#ifndef CORE_BASE_META_HOT_H_
#define CORE_BASE_META_HOT_H_

#include "io/basic_buffer.h"
#include "support/enums.h"
#include "support/magic_constants.h"

namespace tachyon{
namespace core{

struct MetaEntryController {
	typedef MetaEntryController self_type;

	// Ctor
	MetaEntryController(void) :
		gt_available(0),
		gt_anyMissing(0),
		gt_phase(0),
		gt_anyNA(0),
		gt_mixed_phasing(0),
		gt_compression_type(0),
		gt_primtive_type(0),
		biallelic(0),
		simple_snv(0),
		diploid(0),
		mixed_ploidy(0),
		alleles_packed(0),
		all_snv(0)
	{

	}

	MetaEntryController(const self_type& other) :
		gt_available(other.gt_available),
		gt_anyMissing(other.gt_anyMissing),
		gt_phase(other.gt_phase),
		gt_anyNA(other.gt_anyNA),
		gt_mixed_phasing(other.gt_mixed_phasing),
		gt_compression_type(other.gt_compression_type),
		gt_primtive_type(other.gt_primtive_type),
		biallelic(other.biallelic),
		simple_snv(other.simple_snv),
		diploid(other.diploid),
		mixed_ploidy(other.mixed_ploidy),
		alleles_packed(other.alleles_packed),
		all_snv(other.all_snv)
	{}

	void operator=(const self_type& other){
		this->gt_available     = other.gt_available;
		this->gt_anyMissing    = other.gt_anyMissing;
		this->gt_phase         = other.gt_phase;
		this->gt_anyNA         = other.gt_anyNA;
		this->gt_mixed_phasing = other.gt_mixed_phasing;
		this->gt_compression_type = other.gt_compression_type;
		this->gt_primtive_type = other.gt_primtive_type;
		this->biallelic        = other.biallelic;
		this->simple_snv       = other.simple_snv;
		this->diploid          = other.diploid;
		this->mixed_ploidy     = other.mixed_ploidy;
		this->alleles_packed   = other.alleles_packed;
		this->all_snv          = other.all_snv;
	}

	// Dtor
	~MetaEntryController(){}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& entry){
		buffer += (uint16_t)*reinterpret_cast<const uint16_t* const>(&entry);
		return(buffer);

	}

	void operator=(const uint16_t& value){
		const self_type* const other = reinterpret_cast<const self_type* const>(&value);
		this->gt_available     = other->gt_available;
		this->gt_anyMissing    = other->gt_anyMissing;
		this->gt_phase         = other->gt_phase;
		this->gt_anyNA         = other->gt_anyNA;
		this->gt_mixed_phasing = other->gt_mixed_phasing;
		this->gt_compression_type = other->gt_compression_type;
		this->gt_primtive_type = other->gt_primtive_type;
		this->biallelic        = other->biallelic;
		this->simple_snv       = other->simple_snv;
		this->diploid          = other->diploid;
		this->mixed_ploidy     = other->mixed_ploidy;
		this->alleles_packed   = other->alleles_packed;
		this->all_snv          = other->all_snv;
	}

	inline uint16_t toValue(void) const{ return((uint16_t)*reinterpret_cast<const uint16_t* const>(this)); }

	/**< Controller field. The first seven fields describes
	 * genotype-specific information. The remainder bit-fields
	 * describes variant-specific information.
	 * 1) If genotype data is available
	 * 2) If there is are any missing genotypes
	 * 3) The phase of all genotypes if diploid and no mixed phase
	 * 4) If all genotypes share the same phase
	 * 5) What compression method genotypes are stored in
	 * 6) What primitive type the genotypes are encoded as
	 * 7) If this variant is biallelic
	 * 8) If this variant is a simple SNP->SNP
	 * 9) If the genotypes are base diploid
	 * 10) If there is mixed ploidy in the genotypes (i.e. EOV values present)
	 * 11) If allele data is bit-packed
	 * 12) If all alleles are SNVs
	 */
	uint16_t
	    gt_available:        1, // if there is any GT data
        gt_anyMissing:       1, // any missing
        gt_phase:            1, // all phased/unphased
		gt_anyNA:            1, // any NA
        gt_mixed_phasing:    1, // has mixed phasing
		gt_compression_type: 3, // GT compression algorithm
		gt_primtive_type:    2, // type of compression primitive (BYTE, uint16_t, U32, U64)
        biallelic:           1, // is biallelic
        simple_snv:          1, // is simple SNV->SNV
		diploid:             1, // is diploid
		mixed_ploidy:        1, // has mixed ploidy (e.g. X chromosome or CNV)
		alleles_packed:      1, // are the alleles packed in a BYTE
		all_snv:             1; // are all ref/alt simple SNVs
};

}
}

#endif /* CORE_BASE_META_HOT_H_ */
