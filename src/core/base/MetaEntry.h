#ifndef CORE_BASE_METAENTRY_H_
#define CORE_BASE_METAENTRY_H_

#include <limits>

#include "header/Header.h"
#include "MetaHot.h"
#include "MetaCold.h"

namespace Tachyon{
namespace Core{

/**< Envelope record for meta hot-cold split
 * Is used primarily in the simple API. Has a
 * high cost as cold data is copied. Cold
 * meta allele data is always read even though
 * in most cases it is not required
 */
struct MetaEntry{
private:
	typedef MetaHot hot_entry;
	typedef MetaCold cold_entry;
	typedef StreamContainer container_type;
	typedef Header header_type;

public:
	MetaEntry() :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		loaded_cold(false)
	{}

	MetaEntry(const hot_entry& hot) :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		hot(hot),
		loaded_cold(false)
	{}

	MetaEntry(const hot_entry& hot, const cold_entry& cold) :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		hot(hot),
		cold(cold),
		loaded_cold(true)
	{}

	~MetaEntry(){ /* do nothing */ };

	/**<
	 * Translates MetaEntry record into a VCF record
	 * for output to the target ostream
	 * @param dest
	 * @param header
	 * @param blockContigID
	 * @param blockPos
	 */
	void toVCFString(std::ostream& dest, const header_type& header, const S32& blockContigID, const U64& blockPos) const{
		dest.write(&header.getContig(blockContigID).name[0], header.getContig(blockContigID).name.size()) << '\t';
		dest << blockPos + this->hot.position + 1 << '\t';

		// If we have cold meta
		if(this->hasLoadedColdMeta()){
			if(this->cold.n_ID == 0) dest.put('.');
			else dest.write(this->cold.ID, this->cold.n_ID);
			dest << '\t';
			if(this->hot.controller.biallelic && this->hot.controller.simple){
				dest << this->hot.ref_alt.getRef() << '\t' << this->hot.ref_alt.getAlt();
			}
			else {
				dest.write(this->cold.alleles[0].allele, this->cold.alleles[0].l_allele);
				dest << '\t';
				U16 j = 1;
				for(; j < this->cold.n_allele - 1; ++j){
					dest.write(this->cold.alleles[j].allele, this->cold.alleles[j].l_allele);
					dest.put(',');
				}
				dest.write(this->cold.alleles[j].allele, this->cold.alleles[j].l_allele);
			}
			if(std::isnan(this->cold.QUAL)) dest << "\t.\t";
			else dest << '\t' << this->cold.QUAL << '\t';
		}
	}

private:
	inline const bool& hasLoadedColdMeta(void) const{ return(this->loaded_cold); }

public:
	U32 info_pattern_id;
	U32 filter_pattern_id;
	U32 format_pattern_id;

	// Meta objects
	hot_entry hot;
	cold_entry cold;

private:
	bool loaded_cold;
};

}
}

#endif /* CORE_BASE_METAENTRY_H_ */
