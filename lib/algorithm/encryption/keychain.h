#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_H_

#include "containers/variant_block.h"
#include "io/basic_buffer.h"
#include "keychain_key.h"
#include "containers/components/generic_iterator.h"

namespace tachyon{

/**<
 * Primary encryption/decryption keychain for a tachyon
 * archive. The keychain does NOT contain any information
 * to ascertain the relationship between the keychain and
 * the correct archive.
 */
class Keychain {
public:
    typedef Keychain                  self_type;
	typedef std::size_t               size_type;
    typedef KeychainKey*              value_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef containers::VariantBlock  variant_block_type;
    typedef containers::DataContainer container_type;
    typedef std::unordered_map<uint64_t, uint32_t> htable_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    Keychain();
    Keychain(const uint32_t start_capacity);
    Keychain(const self_type& other);
    Keychain(Keychain&& other) = delete;
    Keychain& operator=(Keychain&& other) = delete;
    Keychain& operator=(const self_type& other){
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
    ~Keychain();

    Keychain& operator+=(const self_type& other){
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
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->entries_[0]); }
	inline iterator end(){ return iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries_]); }

	void operator+=(KeychainKey& keychain);

	void resize(const size_type new_capacity);
	void resize(void);

	uint64_t GetRandomHashIdentifier(const bool store = true);
	bool GetHashIdentifier(const uint64_t& value, uint32_t& match);

private:
	bool BuildHashTable(void);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& keychain){
		stream.write(constants::FILE_HEADER.data(), constants::FILE_HEADER_LENGTH);
		stream.write(reinterpret_cast<const char*>(&keychain.n_entries_),  sizeof(size_type));
		stream.write(reinterpret_cast<const char*>(&keychain.n_capacity_), sizeof(size_type));

		// Iteratively write keychain keys.
		for(uint32_t i = 0; i < keychain.size(); ++i)
			stream << *keychain[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& keychain){
		char header[constants::FILE_HEADER_LENGTH];
		stream.read(&header[0], constants::FILE_HEADER_LENGTH);
		if(strncmp(&header[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Illegal keychain..." << std::endl;
			return(stream);
		}

		stream.read(reinterpret_cast<char*>(&keychain.n_entries_),  sizeof(size_type));
		stream.read(reinterpret_cast<char*>(&keychain.n_capacity_), sizeof(size_type));

		for(int i = 0; i < keychain.size(); ++i)
			delete keychain.entries_[i];
		delete [] keychain.entries_;

		keychain.entries_ = new value_type[keychain.n_capacity_];
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

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& keychain){
		buffer += keychain.n_entries_;
		buffer += keychain.n_capacity_;
		for(uint32_t i = 0; i < keychain.size(); ++i) buffer += *keychain[i];
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& keychain){
		// Todo
		/*
		buffer >> keychain.n_entries_;
		buffer >> keychain.n_capacity_;
		delete [] keychain.entries_;
		keychain.entries_ = new value_type[keychain.n_capacity_];
		for(uint32_t i = 0; i < keychain.size(); ++i) buffer >> *keychain.entries_[i];
		keychain.BuildHashTable();
		*/
		return(buffer);
	}

private:
    size_type   n_entries_;
    size_type   n_capacity_;
    pointer     entries_;
    htable_type htable_;
};


}



#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_H_ */
