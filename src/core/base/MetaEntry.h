#ifndef CORE_BASE_METAENTRY_H_
#define CORE_BASE_METAENTRY_H_

#include "MetaHot.h"
#include "MetaCold.h"

namespace Tachyon{
namespace Core{

/**< Envelope record for meta hot-cold split
 * Is used primarily in the simple API. Has a
 * high cost as cold data is copied.
 */
struct MetaEntry{
private:
	typedef MetaHot hot_entry;
	typedef MetaCold cold_entry;
	typedef StreamContainer container_type;

public:
	MetaEntry() : hot(nullptr){}
	MetaEntry(const hot_entry* hot, container_type& container) : hot(hot), cold(&container.buffer_data_uncompressed[hot->virtual_offset_cold_meta]){}
	~MetaEntry(){ /* do nothing */ };

	void operator()(const hot_entry* hot, container_type& container){
		this->hot = hot;
		this->cold(&container.buffer_data_uncompressed[this->hot->virtual_offset_cold_meta]);
	}

public:
	const hot_entry* hot;
	cold_entry cold;
};

}
}

#endif /* CORE_BASE_METAENTRY_H_ */
