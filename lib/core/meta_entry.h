#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "containers/components/variant_block_footer.h"
#include "containers/components/variant_block_header.h"
#include "core/header/variant_header.h"
#include "meta_allele.h"
#include "meta_entry_controller.h"
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
	typedef MetaEntryController         controller_type;

public:
	MetaEntry();
	MetaEntry(const self_type& other);
	MetaEntry(const bcf1_t* record);
	MetaEntry(const bcf1_t* record, const uint64_t position_offset);
	~MetaEntry();

	// Check if a field is set
	inline bool CheckInfoField(const datablock_footer_type& block, const uint32_t info_identifier) const{
		return(block.info_patterns[this->info_pattern_id][info_identifier]);
	}

	inline bool CheckFormatField(const datablock_footer_type& block, const uint32_t format_identifier) const{
		return(block.format_patterns[this->format_pattern_id][format_identifier]);
	}

	inline bool CheckFilterField(const datablock_footer_type& block, const uint32_t filter_identifier) const{
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
	inline const uint16_t& GetNumberAlleles(void){ return(this->n_alleles); }
	inline uint32_t GetContigID(void){ return(this->contigID); }
	inline uint64_t GetPosition(void){ return(this->position); }

	inline const float& GetQuality(void) const{ return(this->quality); }
	inline const std::string& GetName(void) const{ return(this->name); }
	inline const uint16_t& GetNumberAlleles(void) const{ return(this->n_alleles); }
	inline uint32_t GetContigID(void) const{ return(this->contigID); }
	inline uint64_t GetPosition(void) const{ return(this->position); }

	// Set and get for patterns
	inline int32_t& GetInfoPatternId(void){ return(this->info_pattern_id); }
	inline int32_t& GetFormatPatternId(void){ return(this->format_pattern_id); }
	inline int32_t& GetFilterPatternId(void){ return(this->filter_pattern_id); }
	inline const int32_t& GetInfoPatternId(void) const{ return(this->info_pattern_id); }
	inline const int32_t& GetFormatPatternId(void) const{ return(this->format_pattern_id); }
	inline const int32_t& GetFilterPatternId(void) const{ return(this->filter_pattern_id); }

	/**<
	 * Check if it is possible to pack the REF and ALT allele strings into
	 * a single uint8_t. The input allelic data has to be diploid, biallelic and
	 * match the regular expression pattern "^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$"
	 * @return Returns TRUE if it is possible to bitpack data or FALSE otherwise
	 */
	bool UsePackedRefAlt(void) const;

	/**<
	 * Bitpack biallelic, diploid REF and ALT data into a single uint8_t. Failure
	 * to check for validity beforehand with `usePackedRefAlt` may result in
	 * errors.
	 * @return Returns a bitpacked uint8_t
	 */
	uint8_t PackRefAltByte(void) const;

	/**<
	 * Returns the alleles as a concatenated string. For example the alleles
	 * A and T is returned as "A,T". This function is used primarily in the
	 * UpdateHtslibVcfRecord() function for exporting into a htslib bcf1_t record.
	 * @return Returns a concatenated string of alleles.
	 */
	std::string GetAlleleString(void) const{
		std::string ret = this->alleles[0].toString();
		for(uint32_t i = 1; i < this->n_alleles; ++i)
			ret += "," + this->alleles[i].toString();
		return(ret);
	}

	/**<
	 * Updates a htslib bcf1_t record with data available in this meta record.
	 * This function is used when converting yon1_t records to bcf1_t records.
	 * @param rec Input bcf1_t record that has been allocated.
	 * @param hdr Input bcf hdr structure converted from tachyon header.
	 * @return Returns the input bcf1_t record pointer.
	 */
	bcf1_t* UpdateHtslibVcfRecord(bcf1_t* rec, bcf_hdr_t* hdr) const{
		rec->rid = this->contigID;
		rec->pos = this->position;
		bcf_update_id(hdr, rec, this->name.data());
		bcf_update_alleles_str(hdr, rec, this->GetAlleleString().data());
		if(std::isnan(this->quality)) bcf_float_set_missing(rec->qual);
		else rec->qual = this->quality;

		return(rec);
	}

public:
	// Markup: populate from streams
	controller_type controller;
	uint8_t  n_base_ploidy;
	uint16_t   n_alleles;
	int32_t   info_pattern_id;   // Info pattern ID
	int32_t   filter_pattern_id; // Filter pattern ID
	int32_t   format_pattern_id; // Format pattern ID
	float quality;
	uint64_t   contigID;
	uint64_t   position;
	// Todo: Add end_position_longest. This would allow us to directly query
	//       precomputed end positions.
	// Todo: Axtend controller to uint32_t and add fields below or optionally add
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
