#ifndef CORE_BASE_META_HOT_H_
#define CORE_BASE_META_HOT_H_

/*======  Dependencies  ======*/
#include "../io/basic_buffer.h"
#include "../support/enums.h"
#include "../support/MagicConstants.h"

namespace tachyon{
namespace core{

enum TACHYON_GT_PRIMITIVE_TYPE{
	YON_GT_BYTE = 0,
	YON_GT_U16  = 1,
	YON_GT_U32  = 2,
	YON_GT_U64  = 3
};

enum TACHYON_GT_TYPE{
	YON_GT_RLE_DIPLOID_BIALLELIC,
	YON_GT_RLE_DIPLOID_NALLELIC,
	YON_GT_BCF_DIPLOID,
	YON_GT_BCF_STYLE
};

struct MetaHotController{
	typedef MetaHotController self_type;

	// Ctor
	explicit MetaHotController(void) :
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
		gt_available(0),
		alleles_packed(0),
		unused(0)
	{}

	MetaHotController(const self_type& other){
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
		this->gt_available     = other.gt_available;
		this->alleles_packed   = other.alleles_packed;
		this->unused           = other.unused;
	}

	// Dtor
	~MetaHotController(){}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& entry){
		buffer += (U16)*reinterpret_cast<const U16* const>(&entry);
		return(buffer);

	}

	void operator=(const U16& value){
		const self_type* const other = reinterpret_cast<const self_type* const>(&value);
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
		this->gt_available     = other->gt_available;
		this->alleles_packed   = other->alleles_packed;
		this->unused           = other->unused;
	}

	inline const U16 toValue(void) const{ return((U16)*reinterpret_cast<const U16* const>(this)); }

	/**< Controller field. Describes
	 * 1) If any of the genotypes have missing values
	 * 2) If all the genotypes are phased or not
	 * 3) If any genotype has NA/EOV values
	 * 4) If there is mixed phases in the genotypes
	 * 5) If variant site is biallelic
	 * 6) If the genotypes requires "simple" encoding
	 * 7) If the genotypes are run-length encoded or not
	 * 8) What machine word-size the run-length encoded objects are
	 * 9) If the variant site is diploid
	 * 10) If there is mixed ploidy
	 * 11) This variant site have genotypes
	 * 12-16) Reserved for future use
	 *
	 * If the genotypes are NOT diploid then
	 */
	U16 gt_anyMissing:    1, // any missing
        gt_phase:         1, // all phased
		gt_anyNA:         1, // any NA
        gt_mixed_phasing: 1, // has mixed phasing
		gt_compression_type: 3, // uses RLE compression
		gt_primtive_type: 2, // type of RLE (BYTE, U16, U32, U64)
        biallelic:        1, // is biallelic
        simple_snv:       1, // is simple SNV->SNV
		diploid:          1, // is diploid
		mixed_ploidy:     1, // has mixed ploidy (e.g. X chromosome or CNV)
        gt_available:     1, // if there is any GT data
		alleles_packed:   1,
		unused:           1; // reserved
};

}
}

#endif /* CORE_BASE_META_HOT_H_ */
