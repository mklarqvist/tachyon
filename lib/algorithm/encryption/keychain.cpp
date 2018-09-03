#include <openssl/rand.h>

#include "keychain.h"

namespace tachyon{

Keychain::Keychain() :
	n_entries_(0),
	n_capacity_(100000),
	entries_(new value_type[this->n_capacity_])
{
}

Keychain::Keychain(const uint32_t start_capacity) :
	n_entries_(0),
	n_capacity_(start_capacity),
	entries_(new value_type[this->n_capacity_])
{}

Keychain::Keychain(const self_type& other) :
	n_entries_(other.n_entries_),
	n_capacity_(other.n_capacity_),
	entries_(new value_type[this->n_capacity_])
{
	for(uint32_t i = 0; i < this->size(); ++i)
		this->entries_[i] = other.entries_[i]->Clone();
	this->BuildHashTable();
}

Keychain::~Keychain(){
	for(int i = 0; i < this->size(); ++i)
		delete this->entries_[i];
	delete [] this->entries_;
}

void Keychain::operator+=(KeychainKey& keychain){
	if(this->size() + 1 == this->capacity())
		this->resize();

	this->entries_[this->n_entries_++] = keychain.Clone();
}

void Keychain::resize(const size_type new_capacity){
	if(new_capacity < this->capacity()){
		this->n_entries_ = new_capacity;
		return;
	}

	pointer old = this->entries_;
	this->entries_ = new value_type[new_capacity];
	for(uint32_t i = 0; i < this->size(); ++i) this->entries_[i] = old[i];
	delete [] old;
	this->n_capacity_ =  new_capacity;
}

void Keychain::resize(void){ this->resize(this->capacity()*2); }

uint64_t Keychain::GetRandomHashIdentifier(const bool store){
	uint8_t RANDOM_BYTES[32];
	uint64_t value = 0;
	while(true){
		RAND_bytes(&RANDOM_BYTES[0], 32);
		value = XXH64(&RANDOM_BYTES[0], 32, 1337);

		if(value == 0) continue;

		htable_type::const_iterator it = this->htable_.find(value);
		if(it == this->htable_.end()){
			if(store)
				this->htable_[value] = this->size();
			break;
		}
	}

	return value;
}

bool Keychain::GetHashIdentifier(const uint64_t& value, uint32_t& match){
	htable_type::const_iterator it = this->htable_.find(value);
	if(it != this->htable_.end()){
		match = it->second;
		return(true);
	}

	return false;
}

bool Keychain::BuildHashTable(void){
	this->htable_.clear();

	for(uint32_t i = 0; i < this->size(); ++i)
		this->htable_[this->at(i)->field_id] = i;

	return true;
}


}
