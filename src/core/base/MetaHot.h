#ifndef CORE_BASE_METAHOT_H_
#define CORE_BASE_METAHOT_H_

/*======  Dependencies  ======*/
#include "../../io/BasicBuffer.h"
#include "../../support/MagicConstants.h"

namespace Tachyon{
namespace Core{

enum TACHYON_RLE_TYPE{
	YON_BYTE = 0,
	YON_U16 = 1,
	YON_U32 = 2,
	YON_U64 = 3
};

enum TACHYON_GT_TYPE{
	YON_GT_UNKNOWN = -1,
	YON_GT_RLE_DIPLOID_BIALLELIC,
	YON_GT_RLE_DIPLOID_NALLELIC,
	YON_GT_DIPLOID_BCF
};

struct MetaHotController{
	// Ctor
	explicit MetaHotController(void) :
		gt_anyMissing(0),
		gt_phase(0),
		gt_anyNA(0),
		gt_mixed_phasing(0),
		biallelic(0),
		simple(0),
		gt_rle(0),
		gt_primtive_type(0),
		diploid(0),
		mixed_ploidy(0),
		unused(6)
	{}

	// Dtor
	~MetaHotController(){}

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
	 * 11-16) Reserved for future use
	 *
	 * If the genotypes are NOT diploid then
	 */
	U16 gt_anyMissing:    1, // any missing
        gt_phase:         1, // all phased
		gt_anyNA:         1, // any NA
        gt_mixed_phasing: 1, // has mixed phasing
        biallelic:     1, // is biallelic
        simple:        1, // is simple SNV->SNV
        gt_rle:           1, // uses RLE compression
		gt_primtive_type:      2, // type of RLE (BYTE, U16, U32, U64)
		diploid:       1, // is diploid
		mixed_ploidy:  1, // has mixed ploidy (e.g. X chromosome or CNV)
		unused:        5; // reserved
};

/**
 * MetaHotRefAlt:
 * @brief Supportive structure for internal use only. Helper for ref/alt allele encodings
 * Heuristic approach storing the reference/alternative
 * allele information in all cases where the variant site
 * is bi-allelic and simple SNV->SNV change. If this
 * is true then we can store the reference/alternative
 * in a single byte as two nibbles (4 bits). If the variant
 * site does not meet this criterion then the allele data
 * is stored in the cold meta sub-structure.
 */
struct MetaHotRefAlt{
private:
	typedef MetaHotRefAlt self_type;

public:
	MetaHotRefAlt();
	~MetaHotRefAlt();

	inline void operator=(const BYTE& other){
		this->alt = other & 15;
		this->ref = (other >> 4) & 15;
	}

	inline void setMissing(void){
		this->ref = Tachyon::Constants::REF_ALT_N;
		this->alt = Tachyon::Constants::REF_ALT_N;
	}

	inline const char getRef(void) const{ return(Constants::REF_ALT_LOOKUP[this->ref]); }
	inline const char getAlt(void) const{ return(Constants::REF_ALT_LOOKUP[this->alt]); }

	bool setRef(const char& c);
	bool setAlt(const char& c);

public:
	BYTE ref: 4,
	     alt: 4;
};

/**
 * MetaHot:
 * @brief Contains the hot component of the hot-cold split of a
 * variant site meta information
 * Hot sub-structure of a variant sites meta information. This
 * structure requires a CPU that allows non-aligned memory access.
 * Using a packed entry permits the reinterpret_cast of this
 * struct directly from a byte stream.
 *
 */
#pragma pack(push, 1)
struct __attribute__((packed, aligned(1))) MetaHot{
private:
	typedef MetaHot self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef MetaHotController controller_type;
	typedef MetaHotRefAlt allele_type;

public:
	// ctor
	MetaHot();
	MetaHot(const self_type& other);
	MetaHot(self_type&& other) noexcept;
	MetaHot& operator=(const self_type& other) noexcept;
	MetaHot& operator=(self_type&& other) noexcept;
	~MetaHot();

	// Supportive boolean functions
	inline const bool isBiallelic(void) const{ return(this->controller.biallelic); }
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple == true); }
	inline const bool isRLE(void) const{ return(this->controller.gt_rle); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }

	const TACHYON_GT_TYPE getGenotypeType(void) const{
		if(this->controller.gt_rle && this->controller.biallelic && this->controller.diploid) return YON_GT_RLE_DIPLOID_BIALLELIC;
		else if(this->controller.gt_rle && !this->controller.biallelic && this->controller.diploid) return YON_GT_RLE_DIPLOID_NALLELIC;
		else if(!this->controller.gt_rle && this->controller.diploid) return YON_GT_DIPLOID_BCF;
		else return YON_GT_UNKNOWN;
	}

private:
	// Used for debugging only
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' <<
			   (int)*reinterpret_cast<const BYTE* const>(&entry.controller) << '\t' <<
			   entry.ref_alt.getRef() << '\t' << entry.ref_alt.getAlt();
		return(out);
	}

	// Overload operator+= for basic buffer
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		buffer += (U16)*reinterpret_cast<const U16* const>(&entry.controller);
		buffer += (BYTE)*reinterpret_cast<const BYTE* const>(&entry.ref_alt);
		buffer += entry.position;
		return(buffer);
	}

public:
	/**< Controller bit-fields for a variant site */
	controller_type controller;

	/**< Heuristic approach storing the reference/alternative
	 * allele information in all cases where the variant site
	 * is bi-allelic and simple SNV->SNV change. If this
	 * is true then we can store the reference/alternative
	 * in a single byte as two nibbles (4 bits). If the variant
	 * site does not meet this criterion then the allele data
	 * is stored in the cold meta sub-structure.
	 */
	allele_type ref_alt;

	/**< Genomic position in base-0 encoded as the actual
	 * position minus the smallest position in the block
	 * (see BlockEntry.minPosition)
	 */
	U32 position;
};
#pragma pack(pop)

}
}

#endif /* CORE_BASE_METAHOT_H_ */
