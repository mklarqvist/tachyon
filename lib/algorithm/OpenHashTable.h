#ifndef ALGORITHM_OPENHASHTABLE_H_
#define ALGORITHM_OPENHASHTABLE_H_

#include <atomic>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "third_party/xxhash/xxhash.h"
#include "support/type_definitions.h"

namespace tachyon {
namespace hash {

#define HASH_CONSTANT1 452930477

template <class T, class K>
struct OpenHashEntry{
public:
	T key;
	K value;

	OpenHashEntry(){}
};

template <class T, class K>
class HashTable{
public:
    typedef OpenHashEntry<T, K> value_type;

	HashTable(const U64 arraySize);
	HashTable(const U64 arraySize, const U32 n_reprobes);
    ~HashTable();

    // Basic operations
    bool SetItem(const T* key, const K& value, U32 length = sizeof(T));
    bool SetItem(const void* key_address, const T* key, const K& value, U32 length = sizeof(T));
    bool GetItem(const T* key, K*& entry, U32 length = sizeof(T)) const;
    bool GetItem(const void* key_address, const T* key, K*& entry, U32 length = sizeof(T)) const;
    void clear();

    inline const U32& capacity(void) const{return(this->__size);}
    inline const U32& size(void) const{return(this->__occupied);}
    inline bool empty(void) const{return(this->__occupied == 0);}

    inline value_type& operator[](const U32 position){return(*this->__entries[position]);}
    inline const value_type& operator[](const U32 position) const{return(*this->__entries[position]);}
    inline value_type& at(const U32 position){return(*this->__entries[position]);}
    inline const value_type& at(const U32 position) const{return(*this->__entries[position]);}
    inline value_type* pat(const U32 position){return(this->__entries[position]);}
    inline const value_type* pat(const U32 position) const{return(this->__entries[position]);}

protected:
    bool __set(U64& index, const U64& hash2, const T* const key, const K& value);
    bool __get(U64& index, const U64& hash2, const T* const key, K*& value_type) const;

private:
    U64 __occupied;
    U64 __limit;
    U64 __size;
    U32 __retries;
    value_type** __entries;
};

template <class T, class K>
HashTable<T, K>::HashTable(const U64 arraySize) :
	__occupied(0),
	__size(arraySize),
	__retries(50),
	__entries(new value_type*[arraySize])
{
	this->__limit = (U32)((double)arraySize*(1/0.7));
	for(U32 i = 0; i < this->__size; ++i) // Init all to NULL
		this->__entries[i] = nullptr;
}

template <class T, class K>
HashTable<T, K>::HashTable(const U64 arraySize, const U32 n_reprobes) :
	__occupied(0),
	__size(arraySize),
	__retries(n_reprobes),
	__entries(new value_type*[arraySize])
{
	this->__limit = (U32)((double)arraySize*(1/0.7));
	for(U32 i = 0; i < this->__size; ++i) // Init all to NULL
		this->__entries[i] = nullptr;
}

template <class T, class K>
bool HashTable<T, K>::__set(U64& index, const U64& hash2, const T* const key, const K& value){
	U16 RETRIES_COUNTER = 0;
	for(U32 i = 0;;++i){
		if(RETRIES_COUNTER > this->__retries){
			std::cout << "Failed to insert key: " << key << "(" << value << ") after maximum number of reprobes: " << this->__retries << "/" << this->__occupied << std::endl;
			return false;
		}

		index = (index + i*hash2) % this->__size; //double hashing reprobing

		if(this->__entries[index] != nullptr){
			if(this->__entries[index]->key != *key){
				++RETRIES_COUNTER;
				continue;
			} else
				break;

			if(this->__occupied >= this->__limit){ // If we have occupied maximum number of positions in table
				std::cout << "Failed to insert key: " << key << " because the table is full." << std::endl;
				exit(1);
			}
		}
		// Inset
		this->__entries[index]        = new value_type;
		this->__entries[index]->key   = *key;
		this->__entries[index]->value = value;
		this->__occupied++;
		break;
	}
	return true;
}

template <class T, class K>
bool HashTable<T, K>::__get(U64& index, const U64& hash2, const T* const key, K*& entry) const{
	U16 RETRIES_COUNTER = 0;

	for(U32 i = 0;;++i){
		if(RETRIES_COUNTER > this->__retries){ // Guaranteed not be in the table
			return false;
		}
		index = (index + i*hash2) % this->__size; //double hashing reprobing

		if(this->__entries[index] != nullptr){
			if(this->__entries[index]->key != *key){ // If the current position is occupied by other key
				++RETRIES_COUNTER;
				continue;
			} else { // Correct key
				entry = &this->__entries[index]->value;
				return true;
			}
		} else {
			return false; // If we have found a null pointer we know the data does not exist
		}
	}
}


template <class T, class K>
inline bool HashTable<T, K>::SetItem(const T* key, const K &value, U32 length){
	U64 idx = XXH64(key, length, 0);
	const U64 hash2 = XXH64(key, length, HASH_CONSTANT1);
	return(this->__set(idx, hash2, key, value));
}

template <class T, class K>
inline bool HashTable<T, K>::SetItem(const void* key_adress, const T* key, const K &value, U32 length){
	U64 idx = XXH64(key_adress, length, 0);
	const U64 hash2 = XXH64(key_adress, length, HASH_CONSTANT1);
	return(this->__set(idx, hash2, key, value));
}

template <class T, class K>
inline bool HashTable<T, K>::GetItem(const T* key, K*& entry, U32 length) const{
	U64 idx = XXH64(key, length, 0);
	const U64 hash2 = XXH64(key, length, HASH_CONSTANT1);
	return(this->__get(idx, hash2, key, entry));
}

template <class T, class K>
inline bool HashTable<T, K>::GetItem(const void* key_address, const T* key, K*& entry, U32 length) const{
	U64 idx = XXH64(key_address, length, 0);
	const U64 hash2 = XXH64(key_address, length, HASH_CONSTANT1);
	return(this->__get(idx, hash2, key, entry));
}

template <class T, class K>
HashTable<T, K>::~HashTable(){
    for(U32 i = 0; i < this->__size; ++i)
    	delete this->__entries[i];

    delete[] this->__entries;
}

template <class T, class K>
void HashTable<T, K>::clear(void){
	for(U32 i = 0; i < this->__size; ++i){
		delete this->__entries[i];
		this->__entries[i] = nullptr;
	}
	this->__occupied = 0;
}


} /* namespace Hash */
}


#endif /* ALGORITHM_OPENHASHTABLE_H_ */
