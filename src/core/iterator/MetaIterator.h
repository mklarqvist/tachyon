#ifndef CORE_ITERATOR_METAITERATOR_H_
#define CORE_ITERATOR_METAITERATOR_H_

#include "ContainerIterator.h"
#include "../base/MetaEntry.h"
#include "MetaHotIterator.h"
#include "MetaColdIterator.h"

namespace Tachyon{
namespace Iterator{

/**<
 * High-level abstraction of YON meta data. Internally
 * generates MetaEntry records by hard-copy to allow
 * O(1)-time random-access.
 */
class MetaIterator{
private:
	typedef MetaIterator self_type;
	typedef Core::MetaEntry entry_type;
	typedef MetaHotIterator hot_iterator_type;
	typedef Core::MetaHot hot_type;
	typedef MetaColdIterator cold_iterator_type;
	typedef Core::MetaCold cold_type;
	typedef Core::Container container_type;
	typedef ContainerIterator container_iterator_type;

public:
	MetaIterator();
	MetaIterator(const container_type& container);
	MetaIterator(const container_type& container, const Core::TACHYON_GT_TYPE gt_filter_type);
	MetaIterator(const container_type& container_hot, const container_type& container_cold);
	MetaIterator(const container_type& container_hot, const container_type& container_cold, const Core::TACHYON_GT_TYPE gt_filter_type);
	~MetaIterator();

	/**<
	 *
	 * @param gt_type
	 */
	inline void setGenotypeFilter(const Core::TACHYON_GT_TYPE gt_type){ this->gt_filter_option = gt_type; }

	/**<
	 *
	 * @param hot_meta_container
	 * @return
	 */
	bool setup(const container_type& hot_meta_container);

	/**<
	 *
	 * @param hot_meta_container
	 * @param cold_meta_container
	 * @return
	 */
	bool setup(const container_type& hot_meta_container, const container_type& cold_meta_container);

	/**<
	 *
	 * @param info_id_container
	 * @return
	 */
	bool setInfoIDContainer(const container_type& info_id_container);

	/**<
	 *
	 * @param filter_id_container
	 * @return
	 */
	bool setFilterIDContainer(const container_type& filter_id_container);

	/**<
	 *
	 * @param format_id_container
	 * @return
	 */
	bool setFormatIDContainer(const container_type& format_id_container);

	inline entry_type& first(void){ return(this->entries[0]); }

	inline entry_type& last(void){
		if(this->n_entries == 0) return(this->entries[0]);
		return(this->entries[this->n_entries - 1]);
	}

	inline entry_type& operator[](const U32& p){ return(this->entries[p]); }
	inline entry_type& at(const U32& p){ return((*this)[p]); }
	inline entry_type& current(void){ return(this->entries[this->current_position]); }

	inline const S64& size(void) const{ return(this->n_entries); }

	inline void operator++(void){
		if(this->current_position + 1 == this->n_entries) return;
		++this->current_position;
	}

	inline void operator--(void){
		if(this->current_position == 0) return;
		--this->current_position;
	}

	inline void operator+=(const U32& p){ this->current_position += p; }
	inline void operator-=(const U32& p){ this->current_position -= p; }

	/**<
	 * Reset this iterator by reseting all internal
	 * iterators
	 */
	void reset(void);

	inline std::vector<entry_type>& getEntries(void){ return(this->entries); }

private:
	/**< Clears previous data set (if any) */
	inline void clearPrevious(void){ this->entries.clear(); }


private:
	S32  current_position;
	S64  n_entries;
	S64  n_entries_unfiltered;

	Core::TACHYON_GT_TYPE gt_filter_option;

	// Iterators
	hot_iterator_type       hot_iterator;
	cold_iterator_type      cold_iterator;
	container_iterator_type info_id_iterator;
	container_iterator_type filter_id_iterator;
	container_iterator_type format_id_iterator;

	// Entries
	std::vector<entry_type> entries;

	// Container references
	const container_type* container_hot;
	const container_type* container_cold;
	const container_type* info_id_container;
	const container_type* filter_id_container;
	const container_type* format_id_container;
};

}
}

#endif /* CORE_ITERATOR_METAITERATOR_H_ */
