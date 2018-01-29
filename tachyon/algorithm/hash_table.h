#ifndef ALGORITHM_HASH_TABLE_H_
#define ALGORITHM_HASH_TABLE_H_

#include <cstddef>
#include <cstdio>

#include "../support/type_definitions.h"
#include "../third_party/xxhash/xxhash.h"

namespace tachyon{
namespace algorithm{

#define HASH_CONSTANT1 452930477

template <class value_type, class key_type>
struct HashTableEntry{
private:
	typedef HashTableEntry self_type;

public:
	// Rule of 5

public:
	key_type   key;
	value_type value;
};

template <class value_type, class key_type>
class HashTable{
private:
	typedef HashTable                  self_type;
	typedef HashTableEntry<value_type,key_type> value_type;
	typedef value_type&                reference;
	typedef const value_type&          const_reference;
	typedef value_type*                pointer;
	typedef const value_type*          const_pointer;
	typedef std::ptrdiff_t             difference_type;
	typedef std::size_t                size_type;

public:
	HashTable();
	HashTable(const size_t n_capacity);

	// capa
    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const size_type& capacity(void) const{ return(this->n_capacity); }

private:
    size_type n_entries;
    size_type n_capacity;
    pointer   entries_;
};

}
}



#endif /* ALGORITHM_HASH_TABLE_H_ */
