#include <openssl/rand.h>

#include "support/magic_constants.h"
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

Keychain& Keychain::operator=(const self_type& other){
	for(int i = 0; i < this->size(); ++i)
		delete this->entries_[i];
	delete [] this->entries_;

	this->n_entries_ = other.n_entries_;
	this->n_capacity_ = other.n_capacity_;
	this->entries_ = new value_type[this->n_capacity_];

	// Clone data.
	for(uint32_t i = 0; i < this->size(); ++i)
		this->entries_[i] = other.entries_[i]->Clone();

	// Construct hash table for hashed field-id lookups.
	this->BuildHashTable();
	return(*this);
}

Keychain& Keychain::operator+=(const self_type& other){
	if(this->size() + other.size() >= this->capacity()){
		this->resize(this->size() + other.size() + 1024);
	}

	uint32_t offset = 0;
	for(int i = this->size(); i < this->size() + other.size(); ++i, ++offset){
		this->entries_[i] = other.entries_[offset]->Clone();
	}
	this->n_entries_ += other.size();

	// Reconstruct the hash-table for hashed field-id lookups.
	this->BuildHashTable();

	return(*this);
}

void Keychain::operator+=(KeychainKey& keychain){
	this->spinlock_.lock();
	if(this->size() + 1 == this->capacity())
		this->resize();

	this->entries_[this->n_entries_++] = keychain.Clone();
	this->spinlock_.unlock();
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

	// Spinlock for parallel execution without locking with
	// a mutex lock.
	this->spinlock_.lock();
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
	this->spinlock_.unlock();

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

std::ostream& operator<<(std::ostream& stream, const Keychain& keychain){
	stream.write(TACHYON_MAGIC_HEADER.data(), TACHYON_MAGIC_HEADER_LENGTH);
	stream.write(reinterpret_cast<const char*>(&keychain.n_entries_),  sizeof(size_t));
	stream.write(reinterpret_cast<const char*>(&keychain.n_capacity_), sizeof(size_t));

	// Iteratively write keychain keys.
	for(uint32_t i = 0; i < keychain.size(); ++i)
		stream << *keychain[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, Keychain& keychain){
	char header[TACHYON_MAGIC_HEADER_LENGTH];
	stream.read(&header[0], TACHYON_MAGIC_HEADER_LENGTH);
	if(strncmp(&header[0], &TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Illegal keychain..." << std::endl;
		return(stream);
	}

	stream.read(reinterpret_cast<char*>(&keychain.n_entries_),  sizeof(size_t));
	stream.read(reinterpret_cast<char*>(&keychain.n_capacity_), sizeof(size_t));

	for(int i = 0; i < keychain.size(); ++i)
		delete keychain.entries_[i];
	delete [] keychain.entries_;

	keychain.entries_ = new Keychain::value_type[keychain.n_capacity_];
	for(uint32_t i = 0; i < keychain.size(); ++i){
		uint8_t type = 0;
		stream.read(reinterpret_cast<char*>(&type), sizeof(uint8_t));
		if(type == YON_ENCRYPTION_AES_256_GCM)
			keychain[i] = new KeychainKeyGCM<>;
		else {
			std::cerr << "unknown encryption type" << std::endl;
			exit(1);
		}
		stream >> *keychain[i];
	}

	keychain.BuildHashTable();

	return(stream);
}


}
