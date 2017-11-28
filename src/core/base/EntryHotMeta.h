#ifndef CORE_BASE_ENTRYHOTMETA_H_
#define CORE_BASE_ENTRYHOTMETA_H_

/*======  Dependencies  ======*/
#include "../../io/BasicBuffer.h"

namespace Tachyon{
namespace Core{

struct MetaHotController{
	// Ctor
	explicit MetaHotController(void) :
		anyMissing(0),
		allPhased(0),
		mixed_phasing(0),
		biallelic(0),
		simple(0),
		rle(0),
		rle_type(0),
		diploid(0),
		mixed_ploidy(0),
		unused(6)
	{}

	// Dtor
	~MetaHotController(){}

	/**< Controller field. Describes
	 * 1) If any of the genotypes have missing values
	 * 2) If all the genotypes are phased or not
	 * 3) If there is mixed phases in the genotypes
	 * 4) If variant site is biallelic
	 * 5) If the genotypes requires "simple" encoding
	 * 6) If the genotypes are run-length encoded or not
	 * 7) What machine word-size the run-length encoded objects are
	 * 8) If the variant site is diploid
	 * 9) If there is mixed ploidy
	 * 10) Reserved for future use
	 */
	U16 anyMissing:    1, // any missing
        allPhased:     1, // all phased
		mixed_phasing: 1, // has mixed phasing
		biallelic:     1, // is biallelic
		simple:        1, // is simple SNV->SNV
		rle:           1, // uses RLE compression
		rle_type:      2, // type of RLE (BYTE, U16, U32, U64)
		diploid:       1, // is diploid
		mixed_ploidy:  1, // has mixed ploidy (e.g. X chromosome or CNV)
		unused:        6; // reserved
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
	MetaHotRefAlt() : ref(0), alt(0){}
	~MetaHotRefAlt(){}

	void operator=(const BYTE& other){
		this->alt = (other & 15) << 4;
		this->ref |= (other & 15);
	}

	void setMissing(void){
		this->ref = Tachyon::Constants::REF_ALT_N;
		this->alt = Tachyon::Constants::REF_ALT_N;
	}

	bool setRef(const char& c){
		switch(c){
		case 'A': this->ref = Tachyon::Constants::REF_ALT_A; break;
		case 'T': this->ref = Tachyon::Constants::REF_ALT_T; break;
		case 'G': this->ref = Tachyon::Constants::REF_ALT_G; break;
		case 'C': this->ref = Tachyon::Constants::REF_ALT_C; break;
		default:
			std::cerr << Helpers::timestamp("ERROR", "ENCODING") << "Illegal SNV reference..." << std::endl;
			return false;
		}
		return true;
	}

	bool setAlt(const char& c){
		switch(c){
		case 'A': this->alt = Tachyon::Constants::REF_ALT_A; break;
		case 'T': this->alt = Tachyon::Constants::REF_ALT_T; break;
		case 'G': this->alt = Tachyon::Constants::REF_ALT_G; break;
		case 'C': this->alt = Tachyon::Constants::REF_ALT_C; break;
		case 'N': this->alt = Tachyon::Constants::REF_ALT_N; break;
		default:
			std::cerr << Helpers::timestamp("ERROR", "ENCODING") << "Illegal SNV alternative..." << std::endl;
			return false;
		}
		return true;
	}

	inline const char getRef(void) const{ return(Constants::REF_ALT_LOOKUP[this->ref]); }
	inline const char getAlt(void) const{ return(Constants::REF_ALT_LOOKUP[this->alt]); }

	BYTE ref: 4,
	     alt: 4;
};

/**
 * MetaHot:
 * @brief Contains the hot component of the hot-cold split of a variant site meta information
 * Hot sub-structure of a variant sites meta information. This
 * structure requires a CPU that allows non-aligned memory access.
 * Using a packed entry permits the reinterpret_cast of this
 * struct directly from a byte stream.
 *
 */
struct __attribute__((packed)) MetaHot{
private:
	typedef MetaHot self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef MetaHotController controller_type;
	typedef MetaHotRefAlt allele_type;

public:
	// ctor
	MetaHot() :
		position(0),
		AF(0),
		FILTER_map_ID(0),
		INFO_map_ID(0),
		FORMAT_map_ID(0),
		virtual_offset_cold_meta(0),
		virtual_offset_gt(0),
		n_objects(0)
	{}

	// dtor
	~MetaHot(){}

	// Supportive boolean functions
	inline const bool isSingleton(void) const{ return(this->AF == 0); }
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple == true); }
	inline const bool isRLE(void) const{ return(this->controller.rle); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }

	// Supportive functions
	inline const U32& getObjects(void) const{ return(this->n_objects); }

	// Used for debugging only
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' <<
			   (int)*reinterpret_cast<const BYTE* const>(&entry.controller) << '\t' <<
			   entry.ref_alt.getRef() << '\t' << entry.ref_alt.getAlt() << '\t' <<
			   entry.AF << '\t' <<
			   entry.virtual_offset_cold_meta << '\t' <<
			   entry.virtual_offset_gt << '\t' <<
			   entry.n_objects;
		return(out);
	}

	// Overload operator+= for basic buffer
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		buffer += (U16)*reinterpret_cast<const U16* const>(&entry.controller);
		buffer += entry.position;
		buffer += (BYTE)*reinterpret_cast<const BYTE* const>(&entry.ref_alt);
		buffer += entry.AF;
		buffer += entry.FILTER_map_ID;
		buffer += entry.INFO_map_ID;
		buffer += entry.FORMAT_map_ID;
		buffer += entry.virtual_offset_cold_meta;
		buffer += entry.virtual_offset_gt;
		buffer += entry.n_objects;
		return(buffer);
	}

public:
	/**< Controller bit-fields for a variant site */
	controller_type controller;

	/**< Genomic position in base-0 encoded as the actual
	 * position minus the smallest position in the block
	 * (see BlockEntry.minPosition)
	 */
	U32 position;

	/**< Heuristic approach storing the reference/alternative
	 * allele information in all cases where the variant site
	 * is bi-allelic and simple SNV->SNV change. If this
	 * is true then we can store the reference/alternative
	 * in a single byte as two nibbles (4 bits). If the variant
	 * site does not meet this criterion then the allele data
	 * is stored in the cold meta sub-structure.
	 */
	allele_type ref_alt;

	/**< Allele frequency is precomputed as it is frequently
	 * used in several population-genetics approaches
	 */
	float AF;

	/**< Fields describing set-membership to various filter,
	 * info, and format pattern vectors. These patterns are
	 * implicitly encoded in the block header index (see
	 * BlockEntry)
	 */
	U16 FILTER_map_ID;
	U16 INFO_map_ID;
	U16 FORMAT_map_ID;

	/**< This is the virtual pointer offset into the
	 * cold sub-structure of the hot-cold split.
	 */
	U32 virtual_offset_cold_meta;

	/**< Virtual file offset into the appropriate genotype
	 * container stream. What container is targetted is
	 * encoded in the controller.
	 */
	U32 virtual_offset_gt;

	/**< Number of genotype entries encoded */
	U32 n_objects;
};

}
}

#endif /* CORE_BASE_ENTRYHOTMETA_H_ */
