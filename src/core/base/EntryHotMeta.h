#ifndef CORE_BASE_ENTRYHOTMETA_H_
#define CORE_BASE_ENTRYHOTMETA_H_

#include "../../io/BasicBuffer.h"

namespace Tachyon{
namespace Core{

struct MetaHotController{
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
	~MetaHotController(){}

	U16 anyMissing: 1,  // any missing
        allPhased: 1,   // all phased
		mixed_phasing: 1,// has mixed phasing
		biallelic: 1,   // is biallelic
		simple: 1,      // is simple SNV->SNV
		rle: 1,         // uses RLE compression
		rle_type: 2,   // type of RLE (BYTE, U16, U32, U64)
		diploid: 1,    // is diploid
		mixed_ploidy: 1, // has mixed ploidy (e.g. X chromosome or CNV)
		unused: 6; // reserved
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
	typedef MetaHot self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef MetaHotController controller_type;

public:
	MetaHot() :
		position(0),
		ref_alt(0),
		AF(0),
		FILTER_map_ID(0),
		INFO_map_ID(0),
		FORMAT_map_ID(0),
		virtual_offset_cold_meta(0),
		virtual_offset_gt(0),
		n_objects(0)
	{}
	~MetaHot(){}

	inline const bool isSingleton(void) const{ return(this->AF == 0); }
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple == true); }
	inline const bool isRLE(void) const{ return(this->controller.rle); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }
	inline const U32& getRuns(void) const{ return(this->n_objects); }

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' <<
			   (int)*reinterpret_cast<const BYTE* const>(&entry.controller) << '\t' <<
			   (int)entry.ref_alt << '\t' <<
			   entry.AF << '\t' <<
			   entry.virtual_offset_cold_meta << '\t' <<
			   entry.virtual_offset_gt << '\t' <<
			   entry.n_objects;
		return(out);
	}

	// Overload operator+= for basic buffer
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		buffer += *reinterpret_cast<const U16* const>(&entry.controller);
		buffer += entry.position;
		buffer += entry.ref_alt;
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
	U32 position; // is block.minPosition + position

	/**< Heuristic approach storing the reference/alternative
	 * allele information in all cases where the variant site
	 * is bi-allelic and simple SNV->SNV change. If this
	 * is true then we can store the reference/alternative
	 * in a single byte as two nibbles (4 bits). If the variant
	 * site does not meet this criterion then the allele data
	 * is stored in the cold meta sub-structure.
	 */
	BYTE ref_alt;

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
