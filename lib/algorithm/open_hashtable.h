#ifndef OPENHASHTABLE_H_
#define OPENHASHTABLE_H_

#include <atomic>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "third_party/xxhash/xxhash.h"

namespace tachyon {
namespace hash {

#define HASH_CONSTANT1 452930477

template <class key_type, class value_type>
struct OpenHashEntry{
public:
	key_type   key;
	value_type value;

	OpenHashEntry(){}
};

template <class T, class K>
class HashTable{
public:
    typedef OpenHashEntry<T, K> Entry;

	HashTable(const uint64_t arraySize);
    ~HashTable();

    // Basic operations
    void SetItem(const T* key, const K& value, uint32_t length = sizeof(T));
    void SetItem(const void* key_address, const T* key, const K& value, uint32_t length = sizeof(T));
    bool GetItem(const T* key, K*& entry, uint32_t length = sizeof(T));
    bool GetItem(const void* key_address, const T* key, K*& entry, uint32_t length = sizeof(T));
    void clear();
    inline const uint32_t& size(void) const{return this->__size;}
    inline const uint32_t& occupied(void) const{return this->__occupied;}

    Entry& operator[](const uint32_t position){return *this->__entries[position];}
    K& at(const uint32_t position){return this->__entries.at(position);}
    Entry* pat(const uint32_t position){return this->__entries[position];}

protected:
    void __set(uint64_t& index, const uint64_t& hash2, const T* key, const K& value);
    bool __get(uint64_t& index, const uint64_t& hash2, const T* key, K*& entry);

private:
    uint32_t __occupied;
    uint32_t __limit;
    uint32_t __size;
    uint32_t __retries;

    Entry** __entries;
};

template <class T, class K>
HashTable<T, K>::HashTable(const uint64_t arraySize) :
	__occupied(0),
	__size(arraySize),
	__retries(50),
	__entries(new Entry*[arraySize])
{
	this->__limit = (uint32_t)((double)arraySize*(1/0.7));
	for(uint32_t i = 0; i < this->__size; ++i) // Init all to NULL
		this->__entries[i] = nullptr;
}

template <class T, class K>
void HashTable<T, K>::__set(uint64_t& index, const uint64_t& hash2, const T* key, const K& value){
	uint16_t RETRIES_COUNTER = 0;
	for(uint32_t i = 0;;++i){
		if(RETRIES_COUNTER > this->__retries){ //Table is full! Just continue;

//			for(uint32_t i = 0; i < this->__size; ++i)
//				if(this->__entries[idx] != NULL) std::cout << i << "\t" << this->__entries[idx]->key << std::endl;

			exit(1);
		}

		index = (index + i*hash2) % this->__size; //double hashing reprobing
//		std::cout << "set\t\t" << *__key << "\t" << idx << std::endl;

		if(this->__entries[index] != nullptr){
			const T probedKey = this->__entries[index]->key;
			if (probedKey != *key){ // If the current position is occupied by other key
					RETRIES_COUNTER++;
					continue;
			} else break; // If key already exist then simply return

			if(this->__occupied >= this->__limit){ // If we have occupied maximum number of positions in table
				std::cout << "Failed to insert key: " << key << " because the table is full." << std::endl;
				exit(1);
			}
		}
		//std::cout << "new entry" << std::endl;
		this->__entries[index] = new Entry;
		//std::cout << "inserting stuff" << std::endl;

		// Insert into new position
		this->__entries[index]->key = *key;
		this->__entries[index]->value = value;
		this->__occupied++;
		break;
	}
}

template <class T, class K>
bool HashTable<T, K>::__get(uint64_t& index, const uint64_t& hash2, const T* key, K*& entry){
	uint16_t RETRIES_COUNTER = 0;

	for(uint32_t i = 0;;++i){
		if(RETRIES_COUNTER > this->__retries){ //Table is full! Just continue;
			return false;
		}
		index = (index + i*hash2) % this->__size; //double hashing reprobing
//		std::cout << "retrieve\t" << *key << "\t" << index << "\t" << this->__size << std::endl;

		if(this->__entries[index] != nullptr){
//			std::cout << "not null: testing" << std::endl;
			const T probedKey = this->__entries[index]->key;
			if(probedKey != *key){ // If the current position is occupied by other key
				RETRIES_COUNTER++;
//				std::cout << "retrying beacuse" << probedKey << "\t" << key << std::endl;
				continue;
			} else { // Correct key
//				std::cout << "memory location " << &this->__entries[index]->value << std::endl;
				entry = &this->__entries[index]->value;
				return true;
			}

			if(this->__occupied >= this->__limit){ // If we have occupied maximum number of positions in table
//				std::cout << "Failed to insert key: " << key << " because the table is full." << std::endl;
				exit(1);
			}
		} else {
//			std::cout << "not found" << std::endl;
			return false; // If we have found a null pointer we know the data does not exist
		}
	}
}


template <class T, class K>
inline void HashTable<T, K>::SetItem(const T* key, const K &value, uint32_t length){
	uint64_t idx = XXH64(key, length, 0);
	const uint64_t hash2 = XXH64(key, length, HASH_CONSTANT1);
	this->__set(idx, hash2, key, value);
}

template <class T, class K>
inline void HashTable<T, K>::SetItem(const void* key_adress, const T* key, const K &value, uint32_t length){
	uint64_t idx = XXH64(key_adress, length, 0);
	const uint64_t hash2 = XXH64(key_adress, length, HASH_CONSTANT1);
	this->__set(idx, hash2, key, value);
}

template <class T, class K>
inline bool HashTable<T, K>::GetItem(const T* key, K*& entry, uint32_t length){
	uint64_t idx = XXH64(key, length, 0);
	const uint64_t hash2 = XXH64(key, length, HASH_CONSTANT1);
	return this->__get(idx, hash2, key, entry);
}

template <class T, class K>
inline bool HashTable<T, K>::GetItem(const void* key_address, const T* key, K*& entry, uint32_t length){
	uint64_t idx = XXH64(key_address, length, 0);
	const uint64_t hash2 = XXH64(key_address, length, HASH_CONSTANT1);
	return this->__get(idx, hash2, key, entry);
}

template <class T, class K>
HashTable<T, K>::~HashTable(){
    for(uint32_t i = 0; i < this->__size; ++i)
    	if(this->__entries[i] != nullptr)
    		delete this->__entries[i];

    delete[] this->__entries;
}


} /* namespace Hash */
}


#endif /* OPENHASHTABLE_H_ */
