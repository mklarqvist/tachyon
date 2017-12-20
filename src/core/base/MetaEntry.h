#ifndef CORE_BASE_METAENTRY_H_
#define CORE_BASE_METAENTRY_H_

#include <limits>

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

public:
	MetaEntry() :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		hot(nullptr),
		cold(nullptr)
	{}

	MetaEntry(const hot_entry* hot) :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		hot(hot),
		cold(nullptr)
	{}

	MetaEntry(const hot_entry* hot, const cold_entry& cold) :
		info_pattern_id(std::numeric_limits<S32>::min()),
		filter_pattern_id(std::numeric_limits<S32>::min()),
		format_pattern_id(std::numeric_limits<S32>::min()),
		hot(hot),
		cold(cold)
	{}

	~MetaEntry(){ /* do nothing */ };

	inline void operator()(const hot_entry* hot){ this->hot = hot; }
	inline void operator()(const hot_entry* hot, cold_entry& cold){
		this->hot = hot;
		this->cold = cold;
	}

public:
	// Todo: add membership values here
	U32 info_pattern_id;
	U32 filter_pattern_id;
	U32 format_pattern_id;

	// Meta objects
	const hot_entry* hot;
	cold_entry cold;
};

}
}

#endif /* CORE_BASE_METAENTRY_H_ */
