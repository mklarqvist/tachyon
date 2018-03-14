#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "../containers/components/datablock_header.h"
#include "header/variant_header.h"
#include "meta_allele.h"
#include "meta_cold.h"
#include "meta_hot.h"

namespace tachyon{
namespace core{

/**< Envelope record for meta hot-cold split
 * Is used primarily in the simple API. Has a
 * high(ish) cost as cold data has to be copied.
 * Cold meta allele data is always read even
 * though in most cases it is not required
 */
struct MetaEntry{
private:
	typedef MetaHot                     hot_entry;
	typedef MetaCold                    cold_entry;
	typedef containers::DataContainer   container_type;
	typedef containers::DataBlockFooter datablock_footer_type;
	typedef VariantHeader               header_type;
	typedef io::BasicBuffer             buffer_type;
	typedef MetaAllele                  allele_type;
	typedef MetaHotController           controller_type;

public:
	MetaEntry();
	MetaEntry(const hot_entry& hot);
	MetaEntry(const hot_entry& hot, const cold_entry& cold);
	MetaEntry(const hot_entry& hot, char* cold);
	~MetaEntry();

	/**<
	 * Translates MetaEntry record into a VCF string
	 * for output to the target ostream
	 * @param dest
	 * @param header
	 */
	void toVCFString(std::ostream& dest, const header_type& header) const;

	/**<
	 * Translated MetaEntry record into a VCF string
	 * for output to the target buffer
	 * @param dest
	 * @param header
	 */
	void toVCFString(buffer_type& dest, const header_type& header) const;

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
	inline const bool isRLE(void) const{ return(this->controller.gt_rle); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }
	inline const bool isAnyGTMissing(void) const{ return(this->controller.gt_anyMissing); }
	inline const bool isAnyGTSpecial(void) const{ return(this->controller.gt_anyNA); }
	inline const bool isGTMixedPhasing(void) const{ return(this->controller.gt_mixed_phasing); }
	inline const bool getControllerPhase(void) const{ return(this->controller.gt_phase); }

	const TACHYON_GT_TYPE getGenotypeType(void) const{
		if(this->controller.gt_rle && this->controller.biallelic && this->controller.diploid && !this->controller.gt_anyNA) return YON_GT_RLE_DIPLOID_BIALLELIC;
		else if(this->controller.gt_rle && this->controller.diploid) return YON_GT_RLE_DIPLOID_NALLELIC;
		else if(!this->controller.gt_rle && this->controller.diploid) return YON_GT_BCF_DIPLOID;
		else return YON_GT_UNKNOWN;
	}

	const BYTE getPrimitiveWidth(void) const{
		switch(this->controller.gt_primtive_type){
		case(YON_GT_BYTE):  return(1);
		case(YON_GT_U16): return(2);
		case(YON_GT_U32): return(4);
		case(YON_GT_U64): return(8);
		}
		return(0);
	}

	inline const TACHYON_GT_TYPE getGenotypeEncoding(void) const{ return(this->hot.getGenotypeType()); }
	inline const BYTE getGTPrimitiveWidth(void) const{ return(this->hot.getPrimitiveWidth()); }

	inline const float& getQuality(void) const{ return(this->quality); }
	inline const std::string& getName(void) const{ return(this->name); }
	inline const U16& getNumberAlleles(void) const{ return(this->n_alleles); }
	inline const U32 getContigID(void) const{ return(this->contigID); }
	inline const U64 getPosition(void) const{ return(this->position); }

	inline const char getBiallelicAlleleLiteral(const U32 position) const{
		/*
		if(this->isBiallelic()){
			assert(position == 0 || position == 1);
			switch(position){
			case(0): return(this->hot.ref_alt.getRefAlleleLiteral());
			case(1): return(this->hot.ref_alt.getAltAlleleLiteral());
			default: std::cerr << utility::timestamp("ERROR") << "Illegal allele in biallelic!" << std::endl;
			         exit(1);
			}
		} else {
			std::cerr << utility::timestamp("ERROR") << "Variant is not biallelic!" << std::endl;
			exit(1);
		}
		*/
		return('A');
	}

	// Set and get for patterns
	inline S32& getInfoPatternID(void){ return(this->info_pattern_id); }
	inline S32& getFormatPatternID(void){ return(this->format_pattern_id); }
	inline S32& getFilterPatternID(void){ return(this->filter_pattern_id); }
	inline const S32& getInfoPatternID(void) const{ return(this->info_pattern_id); }
	inline const S32& getFormatPatternID(void) const{ return(this->format_pattern_id); }
	inline const S32& getFilterPatternID(void) const{ return(this->filter_pattern_id); }

public:
	bool       loaded_cold;       // Boolean triggered if cold meta object was overloaded
	hot_entry  hot;               // Hot meta object
	cold_entry cold;              // Cold meta object - can be empty

	// Markup: populate from streams
	//BYTE  refalt;
	controller_type controller;
	U32   n_alleles;
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
