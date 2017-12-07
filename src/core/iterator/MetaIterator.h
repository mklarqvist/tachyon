#ifndef CORE_ITERATOR_METAITERATOR_H_
#define CORE_ITERATOR_METAITERATOR_H_

#include "../base/MetaEntry.h"
#include "MetaHotIterator.h"

namespace Tachyon{
namespace Core{
namespace Iterator{

class MetaIterator{
private:
	typedef MetaIterator self_type;
	typedef MetaEntry entry_type;
	typedef MetaHotIterator hot_iterator_type;
	typedef MetaHot hot_type;
	typedef MetaCold cold_type;
	typedef StreamContainer container_type;

	// Function pointers
	typedef const entry_type& (self_type::*operator_function_type)(const U32& p);

public:
	MetaIterator(container_type& container) :
		loadCold(false),
		n_position(0),
		n_entries(0),
		last_modified(-1),
		hot_iterator(container),
		container_hot(&container),
		container_cold(nullptr),
		operator_function(&self_type::__operatorHotOnly)
	{
		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);
	}

	MetaIterator(container_type& container_hot, container_type& container_cold) :
		loadCold(true),
		n_position(0),
		n_entries(0),
		last_modified(-1),
		hot_iterator(container_hot),
		container_hot(&container_hot),
		container_cold(&container_cold),
		operator_function(&self_type::__operatorBoth)
	{
		assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
		this->n_entries = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);
	}


	inline const entry_type& current(void){
		return(this->*operator_function)(this->n_position);
	}

	inline const entry_type& first(void){
		return(this->*operator_function)(0);
	}

	inline const entry_type& last(void){
		if(this->n_entries == 0) return(this->*operator_function)(0);
		return(this->*operator_function)(this->n_entries - 1);
	}

	inline const entry_type& operator[](const U32& p){
		return (this->*operator_function)(p);
	}

	inline const S32& size(void) const{ return(this->n_entries); }

	// No checks are being made
	inline void operator++(void){ ++this->n_position; }
	inline void operator--(void){ --this->n_position; }
	inline void operator+=(const U32& p){ this->n_position += p; }
	inline void operator-=(const U32& p){ this->n_position -= p; }


private:
	inline const entry_type& __operatorHotOnly(const U32 &p){
		if(this->last_modified == p) return(this->entry);
		this->entry(&hot_iterator[p]);
		this->last_modified = p;
		return(this->entry);
	}

	inline const entry_type& __operatorBoth(const U32 &p){
		if(this->last_modified == p) return(this->entry);
		this->entry(&hot_iterator[p], *this->container_cold);
		this->last_modified = p;
		return(this->entry);
	}

public:
	bool loadCold;
	S32 n_position;
	S32 n_entries;
	S32 last_modified;
	entry_type entry;
	hot_iterator_type hot_iterator;
	container_type* container_hot;
	container_type* container_cold;
	operator_function_type operator_function;
};

}
}
}



#endif /* CORE_ITERATOR_METAITERATOR_H_ */
