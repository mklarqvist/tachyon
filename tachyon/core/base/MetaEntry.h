#ifndef CORE_BASE_METAENTRY_H_
#define CORE_BASE_METAENTRY_H_

#include <limits>

#include "header/Header.h"
#include "MetaHot.h"
#include "MetaCold.h"

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
	typedef MetaHot hot_entry;
	typedef MetaCold cold_entry;
	typedef DataContainer container_type;
	typedef Header header_type;
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
	 * @param blockContigID
	 * @param blockPos
	 */
	void toVCFString(std::ostream& dest, const header_type& header, const S32& blockContigID, const U64& blockPos) const;

	/**<
	 * Translated MetaEntry record into a VCF string
	 * for output to the target buffer
	 * @param dest
	 * @param header
	 * @param blockContigID
	 * @param blockPos
	 */
	void toVCFString(buffer_type& dest, const header_type& header, const S32& blockContigID, const U64& blockPos) const;

	// Check if a field is set
	const bool fieldSet(const std::string& field, const Header& header) const{
		core::HeaderMapEntry* entry = nullptr;
		if(header.getEntry("AC", entry)){
			entry->IDX;
		}

		return false;
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

	inline const TACHYON_GT_TYPE getGenotypeEncoding(void) const{ return(this->hot.getGenotypeType()); }
	inline const BYTE getGTPrimitiveWidth(void) const{ return(this->hot.getPrimitiveWidth()); }

	inline const float getQuality(void) const{ return(this->cold.QUAL); }
	inline const std::string getName(void) const{ return(this->cold.getName()); }
	inline const U16& getNumberAlleles(void) const{ return(this->cold.getNumberAlleles()); }

	// Set and get for patterns
	inline U32& getInfoPatternID(void){ return(this->info_pattern_id); }
	inline U32& getFormatPatternID(void){ return(this->format_pattern_id); }
	inline U32& getFilterPatternID(void){ return(this->filter_pattern_id); }
	inline const U32& getInfoPatternID(void) const{ return(this->info_pattern_id); }
	inline const U32& getFormatPatternID(void) const{ return(this->format_pattern_id); }
	inline const U32& getFilterPatternID(void) const{ return(this->filter_pattern_id); }

private:
	inline const bool& hasLoadedColdMeta(void) const{ return(this->loaded_cold); }

private:
	bool loaded_cold;      // Boolean triggered if cold meta object was overloaded
	U32 info_pattern_id;   // Info pattern ID
	U32 filter_pattern_id; // Filter pattern ID
	U32 format_pattern_id; // Format pattern ID
	hot_entry hot;         // Hot meta object
	cold_entry cold;       // Cold meta object - can be empty
};

}
}

#endif /* CORE_BASE_METAENTRY_H_ */
