#ifndef CORE_ITERATOR_METAHOTITERATOR_H_
#define CORE_ITERATOR_METAHOTITERATOR_H_

namespace Tachyon{
namespace Iterator{

class MetaHotIterator{
	typedef MetaHotIterator self_type;
	typedef Core::StreamContainer container_type;
	typedef Core::MetaHot entry_type;

public:
	MetaHotIterator() : n_entries(0), pos(0), entries(nullptr){}
	MetaHotIterator(const container_type& container) : n_entries(0), pos(0), entries(nullptr){
		this->set(container);
	}
	~MetaHotIterator(){}
	bool set(const container_type& container){
		if(container.buffer_data_uncompressed.pointer == 0)
			return false;

		this->clear();
		this->entries = reinterpret_cast<const entry_type* const>(container.buffer_data_uncompressed.data);
		assert((container.buffer_data_uncompressed.pointer % sizeof(entry_type)) == 0);
		this->n_entries = container.buffer_data_uncompressed.pointer / sizeof(entry_type);
		return true;
	}

	inline bool operator()(const container_type& container){
		return(this->set(container));
	}

	void clear(void){
		this->n_entries = 0;
		this->pos = 0;
		entries = nullptr;
	}

	bool nextEntry(const entry_type*& entry){
		if(this->pos == this->n_entries)
			return false;

		entry = &this->entries[this->pos];
		++this->pos;
		return true;
	}

	bool previousEntry(const entry_type*& entry){
		if(this->pos == 0)
			return false;

		entry = &this->entries[this->pos];
		--this->pos;
		return true;
	}

	inline void operator++(void){ if(this->pos != this->n_entries) ++this->pos; }
	inline void operator--(void){ if(this->pos != 0) --this->pos; }
	inline void operator+=(const S32& p){ if(this->pos + p <= this->n_entries) this->pos += p; }
	inline void operator-=(const S32& p){ if(this->pos - p >= 0) this->pos -= p; }

	inline const entry_type& operator[](const U32& p) const{ return(this->entries[p]); }
	inline const size_t size(void) const{ return(this->n_entries); }
	inline const entry_type& first(void) const{ return(this->entries[0]); }
	inline const entry_type& last(void) const{
		if(this->n_entries == 0) return(this->entries[0]);
		return(this->entries[this->n_entries - 1]);
	}

private:
	S32 n_entries;
	S32 pos;
	const entry_type* entries;
};

}
}

#endif /* CORE_ITERATOR_METAHOTITERATOR_H_ */
