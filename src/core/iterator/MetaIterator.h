#ifndef CORE_ITERATOR_METAITERATOR_H_
#define CORE_ITERATOR_METAITERATOR_H_

#include "ContainerIterator.h"
#include "../base/MetaEntry.h"
#include "MetaHotIterator.h"
#include "MetaColdIterator.h"

namespace Tachyon{
namespace Iterator{

/**<
 * High-level abstraction of YON meta data. Has the unfortunate
 * property that it internally generates MetaEntry records by
 * hard-copy.
 */
class MetaIterator{
private:
	typedef MetaIterator self_type;
	typedef Core::MetaEntry entry_type;
	typedef MetaHotIterator hot_iterator_type;
	typedef Core::MetaHot hot_type;
	typedef MetaColdIterator cold_iterator_type;
	typedef Core::MetaCold cold_type;
	typedef Core::StreamContainer container_type;
	typedef ContainerIterator container_iterator_type;

public:
	MetaIterator() :
		current_position(0),
		n_entries(0),
		n_entries_unfiltered(0),
		gt_filter_option(Core::YON_GT_KEEP_ALL),
		//entries(nullptr),
		container_hot(nullptr),
		container_cold(nullptr),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)
	{

	}

	MetaIterator(const container_type& container) :
		current_position(0),
		n_entries(0),
		n_entries_unfiltered(0),
		gt_filter_option(Core::YON_GT_KEEP_ALL),
		hot_iterator(container),
		//entries(nullptr),
		container_hot(&container),
		container_cold(nullptr),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)

	{
		this->setup(container);
	}

	MetaIterator(const container_type& container, const Core::TACHYON_GT_TYPE gt_filter_type) :
		current_position(0),
		n_entries(0),
		n_entries_unfiltered(0),
		gt_filter_option(gt_filter_type),
		hot_iterator(container),
		//entries(nullptr),
		container_hot(&container),
		container_cold(nullptr),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)

	{
		this->setup(container);
	}

	MetaIterator(const container_type& container_hot, const container_type& container_cold) :
		current_position(0),
		n_entries(0),
		n_entries_unfiltered(0),
		gt_filter_option(Core::YON_GT_KEEP_ALL),
		hot_iterator(container_hot),
		cold_iterator(container_cold, hot_iterator.size()),
		//entries(nullptr),
		container_hot(&container_hot),
		container_cold(&container_cold),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)
	{
		this->setup(container_hot, container_cold);
	}

	MetaIterator(const container_type& container_hot, const container_type& container_cold, const Core::TACHYON_GT_TYPE gt_filter_type) :
		current_position(0),
		n_entries(0),
		n_entries_unfiltered(0),
		gt_filter_option(gt_filter_type),
		hot_iterator(container_hot),
		cold_iterator(container_cold, hot_iterator.size()),
		//entries(nullptr),
		container_hot(&container_hot),
		container_cold(&container_cold),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)
	{
		this->setup(container_hot, container_cold);
	}

	~MetaIterator(){
		this->clearPrevious();
	}

	/**<
	 *
	 * @param gt_type
	 */
	inline void setGenotypeFilter(const Core::TACHYON_GT_TYPE gt_type){ this->gt_filter_option = gt_type; }

	/**<
	 *
	 * @param container
	 * @return
	 */
	bool setup(const container_type& hot_meta_container){
		if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
			return false;

		this->container_hot = &hot_meta_container;
		if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
			this->hot_iterator.setup(hot_meta_container);
		else
			this->hot_iterator.setup(hot_meta_container, this->gt_filter_option);

		this->clearPrevious();

		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries = this->hot_iterator.size();
		this->n_entries_unfiltered = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);

		//this->entries = new entry_type*[this->n_entries];
		this->entries.resize(this->n_entries);
		for(U32 i = 0; i < this->n_entries; ++i)
			this->entries[i] = entry_type(this->hot_iterator[i]);

		return true;
	}

	/**<
	 *
	 * @param container_hot
	 * @param container_cold
	 * @return
	 */
	bool setup(const container_type& hot_meta_container, const container_type& cold_meta_container){
		if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
			return false;

		this->container_hot  = &hot_meta_container;
		this->container_cold = &cold_meta_container;
		if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
			this->hot_iterator.setup(hot_meta_container);
		else
			this->hot_iterator.setup(hot_meta_container, this->gt_filter_option);

		this->clearPrevious();

		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries_unfiltered = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);
		this->n_entries = this->hot_iterator.size();
		if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
			this->cold_iterator.setup(cold_meta_container, this->n_entries_unfiltered);
		else
			this->cold_iterator.setup(cold_meta_container, this->n_entries_unfiltered, this->gt_filter_option);

		//this->entries = new entry_type*[this->n_entries];
		this->entries.resize(this->n_entries);
		for(U32 i = 0; i < this->n_entries; ++i)
			this->entries[i] = entry_type(this->hot_iterator[i], this->cold_iterator[i]);

		return true;
	}

	/**<
	 *
	 * @param container
	 * @return
	 */
	bool setInfoIDContainer(const container_type& info_id_container){
		if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
			return false;

		if(this->n_entries <= 0) return false;
		if(this->n_entries_unfiltered <=0) return false;

		this->info_id_container = &info_id_container;
		this->info_id_iterator.setup(info_id_container);
		const void* p;
		if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
			const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
			U32 j = 0;
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->info_id_iterator.getDataIterator()->currentPointer(p);
				if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
					this->entries[j++].getInfoPatternID() = *(S32*)p;
				++this->info_id_iterator;
			}
		} else {
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->info_id_iterator.getDataIterator()->currentPointer(p);
				this->entries[i].getInfoPatternID() = *(S32*)p;
				++this->info_id_iterator;
			}
		}
		return true;
	}

	/**<
	 *
	 * @param container
	 * @return
	 */
	bool setFilterIDContainer(const container_type& filter_id_container){
		if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
			return false;

		if(this->n_entries == 0) return false;
		if(this->n_entries_unfiltered <=0) return false;

		this->filter_id_container = &filter_id_container;
		this->filter_id_iterator.setup(filter_id_container);
		const void* p;
		if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
			const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
			U32 j = 0;
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->filter_id_iterator.getDataIterator()->currentPointer(p);
				if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
					this->entries[j++].getFilterPatternID() = *(S32*)p;
				++this->filter_id_iterator;
			}
		} else {
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->filter_id_iterator.getDataIterator()->currentPointer(p);
				this->entries[i].getFilterPatternID() = *(S32*)p;
				++this->filter_id_iterator;
			}
		}
		return true;
	}

	/**<
	 *
	 * @param container
	 * @return
	 */
	bool setFormatIDContainer(const container_type& format_id_container){
		if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
			return false;

		if(this->n_entries == 0) return false;
		if(this->n_entries_unfiltered <=0) return false;

		this->format_id_container = &format_id_container;
		this->format_id_iterator.setup(format_id_container);
		const void* p;
		if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
			const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
			U32 j = 0;
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->format_id_iterator.getDataIterator()->currentPointer(p);
				if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
					this->entries[j++].getFormatPatternID() = *(S32*)p;
				++this->format_id_iterator;
			}
		} else {
			for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
				this->format_id_iterator.getDataIterator()->currentPointer(p);
				this->entries[i].getFormatPatternID() = *(S32*)p;
				++this->format_id_iterator;
			}
		}
		return true;
	}

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
	void reset(void){
		this->current_position = 0;
		this->hot_iterator.reset();
		this->cold_iterator.reset();
		this->info_id_iterator.reset();
		this->filter_id_iterator.reset();
		this->format_id_iterator.reset();
	}

	inline std::vector<entry_type>& getEntries(void){ return(this->entries); }

private:
	/**< Clears previous data set (if any) */
	inline void clearPrevious(void){
		/*
		for(U32 i = 0; i < this->n_entries; ++i)
			delete this->entries[i];

		delete [] this->entries;
		*/
		this->entries.clear();
	}


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
	//entry_type** entries;
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
