#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_H_

#include "containers/variant_block.h"
#include "io/basic_buffer.h"
#include "keychain_key.h"

namespace tachyon{
namespace encryption{

template < class KeyType = KeychainKeyGCM<> >
class Keychain{
private:
    typedef Keychain                  self_type;
	typedef std::size_t               size_type;
    typedef KeychainKeyGCM<>          value_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef hash::HashTable<U64, U32> hash_table;
    typedef containers::VariantBlock  variant_block_type;
    typedef containers::DataContainer container_type;

public:
    Keychain();
    Keychain(const U32 start_capacity);
    Keychain(const self_type& other);
    ~Keychain();

    class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
	inline reference operator[](const size_type& position){ return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
	inline pointer data(void){ return(this->entries_); }
	inline const_pointer data(void) const{ return(this->entries_); }
	inline reference front(void){ return(this->entries_[0]); }
	inline const_reference front(void) const{ return(this->entries_[0]); }
	inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->entries_[0]); }
	inline iterator end(){ return iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries_]); }

	void operator+=(const value_type& keychain);

	void resize(const size_type new_capacity);
	void resize(void);

	U64 getRandomHashIdentifier();
	bool getHashIdentifier(const U64& value, U32*& match);

private:
	bool __buildHashTable(void);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& keychain){
		stream.write(constants::FILE_HEADER.data(), constants::FILE_HEADER_LENGTH);
		stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MAJOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MINOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_RELEASE), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&keychain.n_entries_),  sizeof(size_type));
		stream.write(reinterpret_cast<const char*>(&keychain.n_capacity_), sizeof(size_type));
		for(U32 i = 0; i < keychain.size(); ++i) stream << keychain[i];
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& keychain){
		char header[constants::FILE_HEADER_LENGTH];
		stream.read(&header[0], constants::FILE_HEADER_LENGTH);
		if(strncmp(&header[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Illegal keychain..." << std::endl;
			return(stream);
		}

		stream.read(reinterpret_cast<char*>(&keychain.version_major_),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&keychain.version_minor_),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&keychain.version_release_), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&keychain.n_entries_),  sizeof(size_type));
		stream.read(reinterpret_cast<char*>(&keychain.n_capacity_), sizeof(size_type));
		delete [] keychain.entries_;
		keychain.entries_ = new value_type[keychain.n_capacity_];
		for(U32 i = 0; i < keychain.size(); ++i) stream >> keychain[i];
		keychain.__buildHashTable();
		return(stream);
	}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& keychain){
		buffer += keychain.n_entries_;
		buffer += keychain.n_capacity_;
		for(U32 i = 0; i < keychain.size(); ++i) buffer += keychain[i];
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& keychain){
		buffer >> keychain.n_entries_;
		buffer >> keychain.n_capacity_;
		delete [] keychain.entries_;
		keychain.entries_ = new value_type[keychain.n_capacity_];
		for(U32 i = 0; i < keychain.size(); ++i) buffer >> keychain.entries_[i];
		keychain.__buildHashTable();
		return(buffer);
	}

private:
	S32 version_major_;
	S32 version_minor_;
	S32 version_release_;
    size_type   n_entries_;
    size_type   n_capacity_;
    pointer     entries_;
    hash_table* htable_identifiers_;
};


template <class KeyType>
Keychain<KeyType>::Keychain() :
	n_entries_(0),
	n_capacity_(100000),
	entries_(new value_type[this->n_capacity_]),
	htable_identifiers_(new hash_table(500000))
{
}

template <class KeyType>
Keychain<KeyType>::Keychain(const U32 start_capacity) :
	n_entries_(0),
	n_capacity_(start_capacity),
	entries_(new value_type[this->n_capacity_]),
	htable_identifiers_(new hash_table(500000))
{}

template <class KeyType>
Keychain<KeyType>::Keychain(const self_type& other) :
	n_entries_(other.n_entries_),
	n_capacity_(other.n_capacity_),
	entries_(new value_type[this->n_capacity_]),
	htable_identifiers_(new hash_table(500000))
{
	for(U32 i = 0; i < this->size(); ++i) this->entries_[i] = other.entries_[i];
	this->__buildHashTable();
}

template <class KeyType>
Keychain<KeyType>::~Keychain(){
	delete [] this->entries_;
	delete this->htable_identifiers_;
}

template <class KeyType>
void Keychain<KeyType>::operator+=(const value_type& keychain){
	if(this->size() + 1 == this->capacity())
		this->resize();

	this->entries_[this->n_entries_++] = keychain;
}

template <class KeyType>
void Keychain<KeyType>::resize(const size_type new_capacity){
	if(new_capacity < this->capacity()){
		this->n_entries_ = new_capacity;
		return;
	}

	pointer old = this->entries_;
	this->entries_ = new value_type[new_capacity];
	for(U32 i = 0; i < this->size(); ++i) this->entries_[i] = old[i];
	delete [] old;
	this->n_capacity_ =  new_capacity;
}

template <class KeyType>
void Keychain<KeyType>::resize(void){ this->resize(this->capacity()*2); }

template <class KeyType>
U64 Keychain<KeyType>::getRandomHashIdentifier(){
	if(this->htable_identifiers_ == nullptr) return false;
	BYTE RANDOM_BYTES[32];
	U64 value = 0;
	U32* match = nullptr;
	while(true){
		RAND_bytes(&RANDOM_BYTES[0], 32);
		value = XXH64(&RANDOM_BYTES[0], 32, 1337);

		if(value == 0) continue;

		if(this->htable_identifiers_->GetItem(&value, match, sizeof(U64)))
			continue;

		if(!this->htable_identifiers_->SetItem(&value, this->size(), sizeof(U64))){
			std::cerr << utility::timestamp("ERROR","HASH") << "Failed to add field identifier..." << std::endl;
			return 0;
		}
		break;
	}

	return value;
}

template <class KeyType>
bool Keychain<KeyType>::getHashIdentifier(const U64& value, U32*& match){
	if(this->htable_identifiers_ == nullptr) return false;
	if(this->htable_identifiers_->GetItem(&value, match, sizeof(U64))){
		return true;
	}
	return false;
}

template <class KeyType>
bool Keychain<KeyType>::__buildHashTable(void){
	delete this->htable_identifiers_;

	if(this->size()*2 < 50e3) this->htable_identifiers_ = new hash_table(50e3);
	else this->htable_identifiers_ = new hash_table(this->size()*2);
	for(U32 i = 0; i < this->size(); ++i){
		if(!this->htable_identifiers_->SetItem(&this->at(i).field_id, i, sizeof(U64))){
			std::cerr << utility::timestamp("ERROR", "KEYCHAIN") << "Failed to generate ID..." << std::endl;
			return false;
		}
	}
	return true;
}


}
}



#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_H_ */
