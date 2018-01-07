#ifndef CORE_ITERATOR_METAITERATOR_H_
#define CORE_ITERATOR_METAITERATOR_H_

#include "ContainerIterator.h"
#include "../base/MetaEntry.h"
#include "MetaHotIterator.h"
#include "MetaColdIterator.h"

namespace Tachyon{
namespace Iterator{

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
		loadCold(false),
		n_position(0),
		n_entries(0),
		entries(nullptr),
		container_hot(nullptr),
		container_cold(nullptr),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)
	{

	}

	MetaIterator(container_type& container) :
		loadCold(false),
		n_position(0),
		n_entries(0),
		hot_iterator(container),
		entries(nullptr),
		container_hot(&container),
		container_cold(nullptr),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)

	{
		this->set(container);
	}

	MetaIterator(container_type& container_hot, container_type& container_cold) :
		loadCold(true),
		n_position(0),
		n_entries(0),
		hot_iterator(container_hot),
		cold_iterator(container_cold, hot_iterator.size()),
		entries(nullptr),
		container_hot(&container_hot),
		container_cold(&container_cold),
		info_id_container(nullptr),
		filter_id_container(nullptr),
		format_id_container(nullptr)
	{
		this->set(container_hot, container_cold);
	}

	~MetaIterator(){
		this->clearPrevious();
	}

	bool set(container_type& container){
		this->container_hot = &container;
		this->hot_iterator.set(container);

		this->clearPrevious();

		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);

		this->entries = new entry_type*[this->n_entries];
		for(U32 i = 0; i < this->n_entries; ++i)
			this->entries[i] = new entry_type(this->hot_iterator[i]);

		return true;
	}

	bool set(container_type& container_hot, container_type& container_cold){
		this->container_hot  = &container_hot;
		this->container_cold = &container_cold;
		this->hot_iterator.set(container_hot);

		this->clearPrevious();

		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);
		this->cold_iterator.set(container_cold, this->n_entries);
		delete [] this->entries;
		this->entries = new entry_type*[this->n_entries];
		for(U32 i = 0; i < this->n_entries; ++i)
			this->entries[i] = new entry_type(this->hot_iterator[i], this->cold_iterator[i]);

		return true;
	}

	inline bool setInfoIDContainer(container_type& container){
		if(this->n_entries == 0) return false;
		this->info_id_container = &container;
		this->info_id_iterator.setup(container);
		const void* p;
		for(U32 i = 0; i < this->n_entries; ++i){
			this->info_id_iterator.getDataIterator().currentPointer(p);
			this->entries[i]->info_pattern_id = *(S32*)p;
			++this->info_id_iterator;
		}
		return true;
	}

	inline bool setFilterIDContainer(container_type& container){
		if(this->n_entries == 0) return false;
		this->filter_id_container = &container;
		this->filter_id_iterator.setup(container);
		const void* p;
		for(U32 i = 0; i < this->n_entries; ++i){
			this->filter_id_iterator.getDataIterator().currentPointer(p);
			this->entries[i]->filter_pattern_id = *(S32*)p;
			++this->filter_id_iterator;
		}
		return true;
	}

	inline bool setFormatIDContainer(container_type& container){
		if(this->n_entries == 0) return false;
		this->format_id_container = &container;
		this->format_id_iterator.setup(container);
		const void* p;
		for(U32 i = 0; i < this->n_entries; ++i){
			this->format_id_iterator.getDataIterator().currentPointer(p);
			this->entries[i]->format_pattern_id = *(S32*)p;
			++this->format_id_iterator;
		}
		return true;
	}

	inline entry_type& first(void){ return(*this->entries[0]); }

	inline entry_type& last(void){
		if(this->n_entries == 0) return(*this->entries[0]);
		return(*this->entries[this->n_entries - 1]);
	}

	inline entry_type& operator[](const U32& p){ return(*this->entries[p]); }

	inline const S32& size(void) const{ return(this->n_entries); }

	// No checks are being made
	inline void operator++(void){ ++this->n_position; }
	inline void operator--(void){ --this->n_position; }
	inline void operator+=(const U32& p){ this->n_position += p; }
	inline void operator-=(const U32& p){ this->n_position -= p; }

	std::vector<entry_type> getEntries(void);

private:
	/**< Clears previous data set (if any) */
	void clearPrevious(void){
		if(this->n_entries){
			for(U32 i = 0; i < this->n_entries; ++i)
				delete this->entries[i];

			delete [] this->entries;
		}
	}


private:
	bool loadCold;
	S32 n_position;
	S32 n_entries;

	// Iterators
	hot_iterator_type hot_iterator;
	cold_iterator_type cold_iterator;
	container_iterator_type info_id_iterator;
	container_iterator_type filter_id_iterator;
	container_iterator_type format_id_iterator;

	// Entries
	entry_type** entries;

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
