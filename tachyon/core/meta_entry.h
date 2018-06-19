#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "../containers/components/variantblock_header.h"
#include "../containers/components/variantblock_footer.h"
#include "header/variant_header.h"
#include "../containers/data_container.h"
#include "meta_allele.h"
#include "../io/bcf/BCFEntry.h"
#include "variant_controller.h"

namespace tachyon{
namespace core{

struct MetaEntry{
public:
	typedef MetaEntry                   self_type;
	typedef containers::DataContainer   container_type;
	typedef containers::VariantBlockFooter datablock_footer_type;
	typedef VariantHeader               header_type;
	typedef io::BasicBuffer             buffer_type;
	typedef MetaAllele                  allele_type;
	typedef VariantController           controller_type;
	typedef bcf::BCFEntry               bcf_entry_type;

public:
	MetaEntry();
	MetaEntry(const self_type& other);
	MetaEntry(const bcf_entry_type& bcf_entry); // transmute data from bcfentry
	MetaEntry(const bcf_entry_type& bcf_entry, const U64 position_offset); // transmute from bcfentry with positional offset
	~MetaEntry();

	// Check if a field is set
	inline const bool check_info_field(const datablock_footer_type& block, const U32 info_identifier) const{
		return(block.info_bit_vectors[this->info_pattern_id][info_identifier]);
	}

	inline const bool check_format_field(const datablock_footer_type& block, const U32 format_identifier) const{
		return(block.format_bit_vectors[this->format_pattern_id][format_identifier]);
	}

	inline const bool check_filter_field(const datablock_footer_type& block, const U32 filter_identifier) const{
		return(block.filter_bit_vectors[this->filter_pattern_id][filter_identifier]);
	}

	// Boolean checks
	// Supportive boolean functions
	inline const bool hasGT(void) const{ return(this->controller.gt_available); }
	inline const bool isBiallelic(void) const{ return(this->controller.biallelic); }
	inline const bool isBiallelicSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple_snv == true); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }
	inline const bool isAnyGTMissing(void) const{ return(this->controller.gt_anyMissing); }
	inline const bool isAnyGTMixedPloidy(void) const{ return(this->controller.gt_anyNA); }
	inline const bool isGTMixedPhasing(void) const{ return(this->controller.gt_mixed_phasing); }
	inline const bool getControllerPhase(void) const{ return(this->controller.gt_phase); }

	inline const TACHYON_GT_ENCODING getGenotypeEncoding(void) const{ return(TACHYON_GT_ENCODING(this->controller.gt_compression_type)); }
	inline const TACHYON_GT_PRIMITIVE_TYPE getGenotypeType(void) const{ return(TACHYON_GT_PRIMITIVE_TYPE(this->controller.gt_primtive_type)); }

	inline const float& getQuality(void){ return(this->quality); }
	inline const std::string& getName(void){ return(this->name); }
	inline const U16& getNumberAlleles(void){ return(this->n_alleles); }
	inline const U32 getContigID(void){ return(this->contigID); }
	inline const U64 getPosition(void){ return(this->position); }

	inline const float& getQuality(void) const{ return(this->quality); }
	inline const std::string& getName(void) const{ return(this->name); }
	inline const U16& getNumberAlleles(void) const{ return(this->n_alleles); }
	inline const U32 getContigID(void) const{ return(this->contigID); }
	inline const U64 getPosition(void) const{ return(this->position); }

	// Set and get for patterns
	inline S32& getInfoPatternID(void){ return(this->info_pattern_id); }
	inline S32& getFormatPatternID(void){ return(this->format_pattern_id); }
	inline S32& getFilterPatternID(void){ return(this->filter_pattern_id); }
	inline const S32& getInfoPatternID(void) const{ return(this->info_pattern_id); }
	inline const S32& getFormatPatternID(void) const{ return(this->format_pattern_id); }
	inline const S32& getFilterPatternID(void) const{ return(this->filter_pattern_id); }

	/**<
	 * Check if it is possible to pack the REF and ALT allele strings into
	 * a single BYTE. The input allelic data has to be diploid, biallelic and
	 * match the regular expression pattern "^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$"
	 * @return Returns TRUE if it is possible to bitpack data or FALSE otherwise
	 */
	const bool usePackedRefAlt(void) const;

	/**<
	 * Bitpack biallelic, diploid REF and ALT data into a single BYTE. Failure
	 * to check for validity beforehand with `usePackedRefAlt` may result in
	 * errors.
	 * @return Returns a bitpacked BYTE
	 */
	const BYTE packRefAltByte(void) const;

public:
	// Markup: populate from streams
	controller_type controller;
	U16   n_alleles;
	S32   info_pattern_id;   // Info pattern ID
	S32   filter_pattern_id; // Filter pattern ID
	S32   format_pattern_id; // Format pattern ID
	float quality;
	U64   contigID;
	U64   position;
	std::string  name;
	allele_type* alleles;
};

}
}

#endif /* CORE_BASE_META_ENTRY_H_ */
