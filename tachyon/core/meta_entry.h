#ifndef CORE_BASE_META_ENTRY_H_
#define CORE_BASE_META_ENTRY_H_

#include <limits>

#include "../containers/core/datablock_header.h"
#include "header/tachyon_header.h"
#include "meta_cold.h"
#include "meta_hot.h"

namespace tachyon{

// Forward declaration for friendship
namespace containers{
	class MetaContainer;
}

namespace core{

/**< Envelope record for meta hot-cold split
 * Is used primarily in the simple API. Has a
 * high(ish) cost as cold data has to be copied.
 * Cold meta allele data is always read even
 * though in most cases it is not required
 */
struct MetaEntry{
private:
	typedef MetaHot hot_entry;
	typedef MetaCold cold_entry;
	typedef containers::DataContainer container_type;
	typedef TachyonHeader header_type;
	typedef io::BasicBuffer buffer_type;

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
	inline const bool check_info_field(const containers::core::DataBlockHeader& block, const U32 info_identifier) const{
		return(block.info_bit_vectors[this->info_pattern_id][info_identifier]);
	}

	inline const bool check_format_field(const containers::core::DataBlockHeader& block, const U32 format_identifier) const{
		return(block.format_bit_vectors[this->format_pattern_id][format_identifier]);
	}

	inline const bool check_filter_field(const containers::core::DataBlockHeader& block, const U32 filter_identifier) const{
		return(block.filter_bit_vectors[this->filter_pattern_id][filter_identifier]);
	}

	// Boolean checks
	inline const bool isBiallelic(void) const{ return(this->hot.isBiallelic()); }
	inline const bool isSimpleSNV(void) const{ return(this->hot.isSimpleSNV()); }
	inline const bool isRLE(void) const{ return(this->hot.isRLE()); }
	inline const bool isDiploid(void) const{ return(this->hot.isDiploid()); }
	inline const bool isMixedPloidy(void) const{ return(this->hot.isMixedPloidy()); }
	inline const bool isAnyGTMissing(void) const{ return(this->hot.isAnyGTMissing()); }
	inline const bool isAnyGTSpecial(void) const{ return(this->hot.isAnyGTSpecial()); }
	inline const bool getControllerPhase(void) const{ return(this->hot.getControllerPhase()); }
	inline const bool isGTMixedPhasing(void) const{ return(this->hot.isGTMixedPhasing()); }
	inline const bool hasGenotypes(void) const{ return(this->hot.hasGT()); }

	inline const TACHYON_GT_TYPE getGenotypeEncoding(void) const{ return(this->hot.getGenotypeType()); }
	inline const BYTE getGTPrimitiveWidth(void) const{ return(this->hot.getPrimitiveWidth()); }

	inline const float getQuality(void) const{ return(this->cold.QUAL); }
	inline const std::string getName(void) const{ return(this->cold.getName()); }
	inline const U16& getNumberAlleles(void) const{ return(this->cold.getNumberAlleles()); }
	inline const U32 getContigID(void) const{ return(this->hot.contigID); }
	inline const U64 getPosition(void) const{ return(this->hot.position); }

	inline const char getBiallelicAlleleLiteral(const U32 position) const{
		if(this->isBiallelic()){
			assert(position == 0 || position == 1);
			switch(position){
			case(0): return(this->hot.ref_alt.getRefAlleleLiteral());
			case(1): return(this->hot.ref_alt.getAltAlleleLiteral());
			default: std::cerr << "Illegal allele in biallelic!" << std::endl; exit(1);
			}
		} else {
			std::cerr << "cannot get this" << std::endl;
			exit(1);
		}
	}

	// Set and get for patterns
	inline S32& getInfoPatternID(void){ return(this->info_pattern_id); }
	inline S32& getFormatPatternID(void){ return(this->format_pattern_id); }
	inline S32& getFilterPatternID(void){ return(this->filter_pattern_id); }
	inline const S32& getInfoPatternID(void) const{ return(this->info_pattern_id); }
	inline const S32& getFormatPatternID(void) const{ return(this->format_pattern_id); }
	inline const S32& getFilterPatternID(void) const{ return(this->filter_pattern_id); }

private:
	inline const bool& hasLoadedColdMeta(void) const{ return(this->loaded_cold); }

public:
	friend containers::MetaContainer;

	bool       loaded_cold;       // Boolean triggered if cold meta object was overloaded
	S32        info_pattern_id;   // Info pattern ID
	S32        filter_pattern_id; // Filter pattern ID
	S32        format_pattern_id; // Format pattern ID
	hot_entry  hot;               // Hot meta object
	cold_entry cold;              // Cold meta object - can be empty
};

}
}

#endif /* CORE_BASE_META_ENTRY_H_ */
