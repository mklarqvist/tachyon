#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "../containers/components/datablock_header.h"
#include "header/variant_header.h"
#include "../containers/datacontainer.h"
#include "meta_allele.h"
#include "../io/bcf/BCFEntry.h"
#include "variant_controller.h"

namespace tachyon{
namespace core{

struct MetaEntry{
public:
	typedef MetaEntry                   self_type;
	typedef containers::DataContainer   container_type;
	typedef containers::DataBlockFooter datablock_footer_type;
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
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple_snv == true); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }
	inline const bool isAnyGTMissing(void) const{ return(this->controller.gt_anyMissing); }
	inline const bool isAnyGTSpecial(void) const{ return(this->controller.gt_anyNA); }
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

	inline const bool isReferenceNONREF(void) const{
		// If the variant site is not biallelic
		if(this->isBiallelic() == false)
			return false;

		int found = 0;
		if((this->alleles[0].l_allele == 9) && (strncmp(this->alleles[0].allele, "<NON_REF>", 9) == 0)){
			++found;
		}

		if((this->alleles[1].l_allele == 9) && (strncmp(this->alleles[1].allele, "<NON_REF>", 9) == 0)){
			++found;
		}

		if(found == 2) return true;
		if(found == 1){
			if(this->alleles[0].l_allele == 1 && this->alleles[1].l_allele == 9) return true;
			if(this->alleles[0].l_allele == 9 && this->alleles[1].l_allele == 1) return true;
			return false;
		}
		return(false);
	}
	/**<
	 *
	 * @return
	 */
	inline const BYTE packedRefAltByte(void) const{
		// <NON_REF>
		BYTE ref_alt = 0;
		if(this->alleles[0].l_allele == 9 && strncmp(this->alleles[0].allele, "<NON_REF>", 9) == 0){
			ref_alt ^= constants::REF_ALT_NON_REF << 4;
		}

		if(this->alleles[1].l_allele == 9 && strncmp(this->alleles[1].allele, "<NON_REF>", 9) == 0){
			ref_alt ^= constants::REF_ALT_NON_REF << 0;
		}

		// Set mock ref-alt if not simple
		if(this->alleles[0].l_allele != 1 || this->alleles[1].l_allele != 1){
			return(255);
		}

		switch(this->alleles[0].allele[0]){
		case 'A': ref_alt ^= constants::REF_ALT_A << 4; break;
		case 'T': ref_alt ^= constants::REF_ALT_T << 4; break;
		case 'G': ref_alt ^= constants::REF_ALT_G << 4; break;
		case 'C': ref_alt ^= constants::REF_ALT_C << 4; break;
		default:
			std::cerr << utility::timestamp("ERROR", "BCF") << "Illegal SNV reference..." << std::endl;
			std::cerr << this->alleles[1].allele << std::endl;
				std::cerr << this->alleles[0].allele << std::endl;
			exit(1);
		}

		switch(this->alleles[1].allele[0]){
		case 'A': ref_alt ^= constants::REF_ALT_A << 0; break;
		case 'T': ref_alt ^= constants::REF_ALT_T << 0; break;
		case 'G': ref_alt ^= constants::REF_ALT_G << 0; break;
		case 'C': ref_alt ^= constants::REF_ALT_C << 0; break;
		case '.': ref_alt ^= constants::REF_ALT_N << 0; break;
		default:
			std::cerr << utility::timestamp("ERROR", "BCF") << "Illegal SNV alt..." << std::endl;
			exit(1);
		}

		return(ref_alt);
	}

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
