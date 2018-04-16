#ifndef ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_
#define ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_

#include <openssl/evp.h>
#include <openssl/aes.h>
#include <openssl/rand.h>
#include "../../containers/variantblock.h"
#include "../../support/enums.h"

namespace tachyon{
namespace encryption{

struct KeychainKey{
public:
	typedef KeychainKey     self_type;
	typedef io::BasicBuffer buffer_type;

public:
	KeychainKey() : fieldIdentifier(0), encryption_type(YON_ENCRYPTION_NONE){}
	virtual ~KeychainKey(){}

	KeychainKey(const KeychainKey& other) :
		fieldIdentifier(other.fieldIdentifier),
		encryption_type(other.encryption_type)
	{
		memcpy(&this->key[0], &other.key[0], 32);
		memcpy(&this->iv[0],  &other.iv[0],  16);
	}

	KeychainKey& operator=(const KeychainKey& other){
		this->fieldIdentifier = other.fieldIdentifier;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], 32);
		memcpy(&this->iv[0],  &other.iv[0],  16);
		return(*this);
	}

	virtual void print(void){
		for(U32 i = 0; i < 32; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(U32 i = 0; i < 16; ++i) std::cerr << std::hex << (int)this->iv[i];
	}

	inline bool initiateRandom(void){
		RAND_bytes(&this->key[0], 32);
		RAND_bytes(&this->iv[0],  16);
		return true;
	}

private:
	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.fieldIdentifier;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], 32);
		buffer.Add((const char*)&key.iv[0],  16);
		return(buffer);
	}

public:
	U64  fieldIdentifier;
	BYTE encryption_type;
	BYTE key[32];  // 256 bit key
	BYTE iv[16];   // 128 bit initiation vector
};

struct KeychainKeyGCM : public KeychainKey{
public:
	typedef KeychainKeyGCM self_type;
	typedef KeychainKey    parent_type;

public:
	KeychainKeyGCM(){}
	~KeychainKeyGCM(){}

	KeychainKeyGCM(const KeychainKeyGCM& other) :
		parent_type(other)
	{
		memcpy(&this->tag[0], &other.tag[0], 16);
	}

	KeychainKeyGCM& operator=(const KeychainKeyGCM& other){
		this->fieldIdentifier = other.fieldIdentifier;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], 32);
		memcpy(&this->iv[0],  &other.iv[0],  16);
		memcpy(&this->tag[0], &other.tag[0], 16);
		return(*this);
	}

	void print(void){
		for(U32 i = 0; i < 32; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(U32 i = 0; i < 16; ++i) std::cerr << std::hex << (int)this->iv[i];
		std::cerr << '\t';
		for(U32 i = 0; i < 16; ++i) std::cerr << std::hex << (int)this->tag[i];
		std::cerr << std::dec;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& key){
		stream.write(reinterpret_cast<const char*>(&key.fieldIdentifier), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&key.encryption_type), sizeof(BYTE));
		stream.write((char*)&key.key[0], 32);
		stream.write((char*)&key.iv[0],  16);
		stream.write((char*)&key.tag[0], 16);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& key){
		stream.read(reinterpret_cast<char*>(&key.fieldIdentifier), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(BYTE));
		stream.read((char*)&key.key[0], 32);
		stream.read((char*)&key.iv[0],  16);
		stream.read((char*)&key.tag[0], 16);
		return(stream);
	}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.fieldIdentifier;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], 32);
		buffer.Add((const char*)&key.iv[0],  16);
		buffer.Add((const char*)&key.tag[0], 16);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& key){
		buffer >> key.fieldIdentifier;
		buffer >> key.encryption_type;
		buffer.read((char*)&key.key[0], 32);
		buffer.read((char*)&key.iv[0],  16);
		buffer.read((char*)&key.tag[0], 16);
		return(buffer);
	}

public:
	BYTE tag[16];  // validation tag for gcm
};

class Keychain{
private:
    typedef Keychain           self_type;
	typedef std::size_t        size_type;
    typedef KeychainKeyGCM     value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef hash::HashTable<U64, U32> hash_table;
    typedef containers::VariantBlock  variant_block_type;
    typedef containers::DataContainer container_type;

public:
    Keychain() :
    	n_entries_(0),
		n_capacity_(100000),
		entries_(new value_type[this->n_capacity_]),
		htable_identifiers_(new hash_table(500000))
	{
    }

    Keychain(const U32 start_capacity) :
    	n_entries_(0),
		n_capacity_(start_capacity),
		entries_(new value_type[this->n_capacity_]),
		htable_identifiers_(new hash_table(500000))
    {}

    Keychain(const self_type& other) :
    	n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		entries_(new value_type[this->n_capacity_]),
		htable_identifiers_(new hash_table(500000))
    {
    	for(U32 i = 0; i < this->size(); ++i) this->entries_[i] = other.entries_[i];
    	this->__buildHashTable();
    }

    ~Keychain(){
    	delete [] this->entries_;
    	delete this->htable_identifiers_;
    }

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

	inline void operator+=(const value_type& keychain){
		if(this->size() + 1 == this->capacity())
			this->resize();

		this->entries_[this->n_entries_++] = keychain;
	}

	void resize(const size_type new_capacity){
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

	void resize(void){ this->resize(this->capacity()*2); }


	U64 getRandomHashIdentifier(){
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

	bool getHashIdentifier(const U64& value, U32*& match){
		if(this->htable_identifiers_ == nullptr) return false;
		if(this->htable_identifiers_->GetItem(&value, match, sizeof(U64))){
			return true;
		}
		return false;
	}

private:
	bool __buildHashTable(void){
		delete this->htable_identifiers_;

		if(this->size()*2 < 50e3) this->htable_identifiers_ = new hash_table(50e3);
		else this->htable_identifiers_ = new hash_table(this->size()*2);
		for(U32 i = 0; i < this->size(); ++i){
			if(!this->htable_identifiers_->SetItem(&this->at(i).fieldIdentifier, i, sizeof(U64))){
				std::cerr << utility::timestamp("ERROR", "KEYCHAIN") << "Failed to generate ID..." << std::endl;
				return false;
			}
		}
		return true;
	}

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
			std::cerr << "ILLEGAL KEY FILE" << std::endl;
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

class EncryptionDecorator{
public:
	typedef EncryptionDecorator       self_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer stream_container;
	typedef io::BasicBuffer           buffer_type;
	typedef KeychainKeyGCM            aes256gcm_type;
	typedef Keychain                  keychain_type;

public:
	bool encrypt(variant_block_type& block, keychain_type& keychain, TACHYON_ENCRYPTION encryption_type){
		if(encryption_type == YON_ENCRYPTION_AES_256_GCM){
			return(this->encryptAES256(block, keychain));
		}
		return(false);
	}

	bool decryptAES256(variant_block_type& block, keychain_type& keychain){
		if(!this->decryptAES256(block.meta_contig_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_positions_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_refalt_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_controller_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_quality_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_names_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_rle8_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_rle16_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_rle32_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_rle64_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_alleles_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_simple8_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_simple16_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_simple32_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_simple64_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.gt_support_data_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_info_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_filter_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		if(!this->decryptAES256(block.meta_format_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(!this->decryptAES256(block.info_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		}

		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
			if(!this->decryptAES256(block.format_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
		}

		return(true);
	}

	bool encryptAES256(variant_block_type& block, keychain_type& keychain){
		BYTE RANDOM_BYTES[32];
		RAND_bytes(&RANDOM_BYTES[0], 32);
		block.header.blockID = XXH64(&RANDOM_BYTES[0], 32, 1337);

		if(!this->encryptAES256(block.meta_contig_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_positions_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_refalt_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_controller_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_quality_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_names_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_rle8_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_rle16_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_rle32_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_rle64_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_alleles_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_simple8_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_simple16_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_simple32_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_simple64_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.gt_support_data_container, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_info_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_filter_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		if(!this->encryptAES256(block.meta_format_map_ids, keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(!this->encryptAES256(block.info_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		}

		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
			if(!this->encryptAES256(block.format_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		}

		return(true);
	}

	bool encryptAES256(stream_container& container, keychain_type& keychain){
		KeychainKeyGCM entry;
		entry.initiateRandom();
		entry.encryption_type = YON_ENCRYPTION_AES_256_GCM;

		EVP_CIPHER_CTX *ctx = NULL;
		S32 len = 0;

		if(!(ctx = EVP_CIPHER_CTX_new())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
			return false;
		}

		if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption..." << std::endl;
			return false;
		}

		// 16 bytes = 128 bits
		if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption tag..." << std::endl;
			return false;
		}

		if(1 != EVP_EncryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption key and IV..." << std::endl;
			return false;
		}

		// Mask header in encrypted message
		this->buffer.reset();
		io::BasicBuffer temp(65536);
		temp << container.header;
		this->buffer.resize(container.buffer_data.size() + container.buffer_strides.size() + temp.size() + 65536);

		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)this->buffer.data(), &len, (BYTE*)temp.data(), temp.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}
		this->buffer.n_chars = len;

		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)&this->buffer[this->buffer.n_chars], &len, (BYTE*)container.buffer_data.data(), container.buffer_data.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)&this->buffer[this->buffer.n_chars], &len, (BYTE*)container.buffer_strides.data(), container.buffer_strides.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		if(1 != EVP_EncryptFinal_ex(ctx, (BYTE*)this->buffer.data() + len, &len)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to finalise the encryption..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		container.buffer_data.resize(this->buffer.size());
		if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, &entry.tag[0])){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to retrieve the GCM tag..." << std::endl;
			return false;
		}

		EVP_CIPHER_CTX_free(ctx);

		// Trigger encryption flag
		container.header.reset(); // reset data
		container.header.data_header.controller.encryption = YON_ENCRYPTION_AES_256_GCM;
		memcpy(container.buffer_data.data(), this->buffer.data(), this->buffer.size());
		container.buffer_data.n_chars = this->buffer.size();
		container.header.data_header.eLength = this->buffer.size();

		const U64 hashID = keychain.getRandomHashIdentifier();
		container.header.identifier = hashID;
		entry.fieldIdentifier = hashID;
		keychain += entry; // add key to keychain

		return(true);
	}

	bool decryptAES256(stream_container& container, keychain_type& keychain){
		if(container.buffer_data.size() == 0)
			return true;

		if(container.header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
			return true;

		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_AES_256_GCM){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Wrong decryption algorithm used..." << std::endl;
			return false;
		}

		U32* match = nullptr;
		if(keychain.getHashIdentifier(container.header.identifier, match) == false){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Did not find ID in keychain..." << std::endl;
			return false;
		}
		KeychainKeyGCM& entry = keychain[*match];

		EVP_CIPHER_CTX *ctx = NULL;
		S32 len = 0, plaintext_len = 0, ret;

		/* Create and initialise the context */
		if(!(ctx = EVP_CIPHER_CTX_new())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
			return(false);
		}

		/* Initialise the decryption operation. */
		if(!EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the decryption..." << std::endl;
			return(false);
		}

		if(!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the IV context..." << std::endl;
			return(false);
		}

		/* Initialise key and IV */
		if(!EVP_DecryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the decryption..." << std::endl;
			return(false);
		}

		this->buffer.reset();
		if(container.buffer_data.size()){
			this->buffer.resize(container.buffer_data.size() + 65536);
			if(!EVP_DecryptUpdate(ctx, (BYTE*)this->buffer.data(), &len, (BYTE*)container.buffer_data.data(), container.buffer_data.size())){
				std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the decryption..." << std::endl;
				return(false);
			}
			plaintext_len = len;
		}

		/* Set expected tag value. Works in OpenSSL 1.0.1d and later */
		if(!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, entry.tag)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to set expected TAG value..." << std::endl;
			return(false);
		}

		/* Finalise the decryption. A positive return value indicates success,
		 * anything else is a failure - the plaintext is not trustworthy.
		 */
		ret = EVP_DecryptFinal_ex(ctx, (BYTE*)this->buffer.data() + len, &len);
		EVP_CIPHER_CTX_free(ctx);

		if(ret > 0){
			plaintext_len += len;
			this->buffer >> container.header; // unmask encrypted header

			container.buffer_data.reset();
			container.buffer_strides.reset();
			container.buffer_data.resize(container.header.data_header.cLength + 65536);
			container.buffer_strides.resize(container.header.stride_header.cLength + 65536);
			memcpy(container.buffer_data.data(), &this->buffer[this->buffer.iterator_position_], container.header.data_header.cLength);
			memcpy(container.buffer_strides.data(), &this->buffer[this->buffer.iterator_position_ + container.header.data_header.cLength], container.header.stride_header.cLength);
			container.header.data_header.controller.encryption = YON_ENCRYPTION_NONE;
			container.buffer_data.n_chars = container.header.data_header.cLength;
			container.buffer_strides.n_chars = container.header.stride_header.cLength;
			container.header.data_header.eLength = 0;
			return(true);
		} else {
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to validate decryption..." << std::endl;
			return(false);
		}
	}


public:
	buffer_type buffer;
};

}
}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_ */
