#ifndef CORE_ITERATOR_METACOLDITERATOR_H_
#define CORE_ITERATOR_METACOLDITERATOR_H_

#include <cassert>

namespace Tachyon{
namespace Core{
namespace Iterator{

/** Iterator class for MetaCold
 *  Unlike many other iterators in Tachyon this one
 *  has to do a hard copy of the underlying data
 *  from the passed byte stream.
 */
class MetaColdIterator{
	typedef MetaColdIterator self_type;
	typedef Core::StreamContainer container_type;
	typedef Core::MetaCold entry_type;

public:
	MetaColdIterator() : n_entries(0), pos(0), offsets(nullptr), container(nullptr){}

	MetaColdIterator(container_type& container, const U32 n_entries) :
		n_entries(n_entries),
		pos(0),
		offsets(nullptr),
		container(nullptr)
	{
		this->set(container, n_entries);
	}

	~MetaColdIterator(){
		// entries are copied from the data stream
		delete [] this->offsets;
	}

	/**<
	 * Interprets and copies data from the passed container
	 * structure into an internal array of structures
	 * @param  container Data container
	 * @param  n_entries Number of entries in 'container'. This is always provided from the 'MetaHotIterator'
	 * @return Returns TRUE if there is data or FALSE if there is not
	 */
	bool set(container_type& container, const S32 n_entries){
		if(container.buffer_data_uncompressed.pointer == 0)
			return false;

		this->clear();

		delete [] this->offsets;
		this->offsets = new U32[n_entries];
		this->n_entries = n_entries;
		this->container = &container;

		// Todo: we should just calculate the offset positions as this point!
		// Then in the MetaEntry struct we invoke the overloaded () operator in
		// MetaCold
		U32 pos = 0;
		this->offsets[0] = 0;
		pos = *reinterpret_cast<const U32* const>(&container.buffer_data_uncompressed.data[0]);
		for(S32 i = 1; i < this->n_entries; ++i){
			const U32& l_body = *reinterpret_cast<const U32* const>(&container.buffer_data_uncompressed.data[pos]);
			this->offsets[i] = pos;
			//std::cerr << i << '/' << this->n_entries << std::endl;
			pos += l_body;
		}

		assert(pos == container.buffer_data_uncompressed.pointer);
		return(true);
	}

	inline bool operator()(container_type& container, const U32 n_entries){
		return(this->set(container, n_entries));
	}

	void clear(void){
		this->n_entries = 0;
		this->pos = 0;
		delete [] this->offsets;
		this->offsets = nullptr;
		this->container = nullptr;
	}

	bool nextEntry(entry_type& entry){
		if(this->pos == this->n_entries)
			return false;

		this->internal_entry(&this->container->buffer_data_uncompressed.data[this->offsets[this->pos]]);
		entry = this->internal_entry;
		++this->pos;
		return true;
	}

	bool previousEntry(entry_type& entry){
		if(this->pos == 0)
			return false;

		--this->pos;
		this->internal_entry(&this->container->buffer_data_uncompressed.data[this->offsets[this->pos]]);
		entry = this->internal_entry;
		return true;
	}

	inline entry_type& currentEntry(entry_type& entry){ return(this->current(entry)); }
	inline entry_type& current(entry_type& entry){ return(this->internal_entry); }

	inline void operator++(void){ if(this->pos != this->n_entries) ++this->pos; }
	inline void operator--(void){ if(this->pos != 0) --this->pos; }
	inline void operator+=(const S32& p){ if(this->pos + p <= this->n_entries) this->pos += p; }
	inline void operator-=(const S32& p){ if(this->pos - p >= 0) this->pos -= p; }

	inline entry_type operator[](const U32& p){ return(entry_type(&this->container->buffer_data_uncompressed.data[this->offsets[p]])); }
	inline const size_t size(void) const{ return(this->n_entries); }
	inline entry_type first(void){ return(entry_type(&this->container->buffer_data_uncompressed.data[this->offsets[0]])); }
	inline entry_type last(void){
		if(this->n_entries == 0) return(this->first());
		return(entry_type(&this->container->buffer_data_uncompressed.data[this->offsets[this->n_entries - 1]]));
	}

private:
	S32 n_entries;
	S32 pos;
	entry_type internal_entry;
	U32* offsets;
	container_type* container;
};

}
}
}


#endif /* CORE_ITERATOR_METACOLDITERATOR_H_ */
