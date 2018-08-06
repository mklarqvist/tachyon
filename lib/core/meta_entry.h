#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "containers/components/variant_block_footer.h"
#include "containers/components/variant_block_header.h"
#include "core/header/variant_header.h"
#include "meta_allele.h"
#include "variant_controller.h"
#include "containers/data_container.h"

#include "htslib/vcf.h"

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

public:
	MetaEntry();
	MetaEntry(const self_type& other);
	MetaEntry(const bcf1_t* record);
	MetaEntry(const bcf1_t* record, const U64 position_offset);
	~MetaEntry();

	// Check if a field is set
	inline bool CheckInfoField(const datablock_footer_type& block, const U32 info_identifier) const{
		return(block.info_patterns[this->info_pattern_id][info_identifier]);
	}

	inline bool CheckFormatField(const datablock_footer_type& block, const U32 format_identifier) const{
		return(block.format_patterns[this->format_pattern_id][format_identifier]);
	}

	inline bool CheckFilterField(const datablock_footer_type& block, const U32 filter_identifier) const{
		return(block.filter_patterns[this->filter_pattern_id][filter_identifier]);
	}

	// Boolean checks
	// Supportive boolean functions
	inline bool HasGT(void) const{ return(this->controller.gt_available); }
	inline bool IsBiallelic(void) const{ return(this->controller.biallelic); }
	inline bool IsBiallelicSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple_snv == true); }
	inline bool IsDiploid(void) const{ return(this->controller.diploid); }
	inline bool IsMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }
	inline bool IsAnyGTMissing(void) const{ return(this->controller.gt_anyMissing); }
	inline bool IsAnyGTMixedPloidy(void) const{ return(this->controller.gt_anyNA); }
	inline bool IsGTMixedPhasing(void) const{ return(this->controller.gt_mixed_phasing); }
	inline bool GetControllerPhase(void) const{ return(this->controller.gt_phase); }

	inline TACHYON_GT_ENCODING GetGenotypeEncoding(void) const{ return(TACHYON_GT_ENCODING(this->controller.gt_compression_type)); }
	inline TACHYON_GT_PRIMITIVE_TYPE GetGenotypeType(void) const{ return(TACHYON_GT_PRIMITIVE_TYPE(this->controller.gt_primtive_type)); }

	inline const float& GetQuality(void){ return(this->quality); }
	inline const std::string& GetName(void){ return(this->name); }
	inline const U16& GetNumberAlleles(void){ return(this->n_alleles); }
	inline U32 GetContigID(void){ return(this->contigID); }
	inline U64 GetPosition(void){ return(this->position); }

	inline const float& GetQuality(void) const{ return(this->quality); }
	inline const std::string& GetName(void) const{ return(this->name); }
	inline const U16& GetNumberAlleles(void) const{ return(this->n_alleles); }
	inline U32 GetContigID(void) const{ return(this->contigID); }
	inline U64 GetPosition(void) const{ return(this->position); }

	// Set and get for patterns
	inline S32& GetInfoPatternId(void){ return(this->info_pattern_id); }
	inline S32& GetFormatPatternId(void){ return(this->format_pattern_id); }
	inline S32& GetFilterPatternId(void){ return(this->filter_pattern_id); }
	inline const S32& GetInfoPatternId(void) const{ return(this->info_pattern_id); }
	inline const S32& GetFormatPatternId(void) const{ return(this->format_pattern_id); }
	inline const S32& GetFilterPatternId(void) const{ return(this->filter_pattern_id); }

	/**<
	 * Check if it is possible to pack the REF and ALT allele strings into
	 * a single BYTE. The input allelic data has to be diploid, biallelic and
	 * match the regular expression pattern "^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$"
	 * @return Returns TRUE if it is possible to bitpack data or FALSE otherwise
	 */
	bool UsePackedRefAlt(void) const;

	/**<
	 * Bitpack biallelic, diploid REF and ALT data into a single BYTE. Failure
	 * to check for validity beforehand with `usePackedRefAlt` may result in
	 * errors.
	 * @return Returns a bitpacked BYTE
	 */
	BYTE PackRefAltByte(void) const;

public:
	// Markup: populate from streams
	controller_type controller;
	BYTE  n_base_ploidy;
	U16   n_alleles;
	S32   info_pattern_id;   // Info pattern ID
	S32   filter_pattern_id; // Filter pattern ID
	S32   format_pattern_id; // Format pattern ID
	float quality;
	U64   contigID;
	U64   position;
	// Todo: Add end_position_longest. This would allow us to directly query
	//       precomputed end positions.
	// Todo: Axtend controller to U32 and add fields below or optionally add
	//       another controller field with variant-identifying information.
	//       This would support queries on type directly from precomputed
	//       lookups.
	// Fields of interest:
	// isSNV
	// isComplex
	// isIndel
	// isSV
	std::string  name;
	allele_type* alleles;
};

}
}

#endif /* CORE_BASE_META_ENTRY_H_ */
