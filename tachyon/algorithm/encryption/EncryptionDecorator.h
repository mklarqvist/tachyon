#ifndef ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_
#define ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_

#include <openssl/evp.h>
#include <openssl/aes.h>
#include <openssl/rand.h>
#include "../../containers/variantblock.h"
#include "../../support/enums.h"

namespace tachyon{
namespace encryption{

struct KeychainEntry{
public:
	typedef KeychainEntry   self_type;
	typedef io::BasicBuffer buffer_type;

public:
	KeychainEntry(){}
	virtual ~KeychainEntry(){}

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

public:
	BYTE key[32];  // 256 bit key
	BYTE iv[16];   // 128 bit initiation vector
};

struct KeychainEntryGCM : public KeychainEntry{
public:
	typedef KeychainEntryGCM self_type;

public:
	KeychainEntryGCM(){}
	~KeychainEntryGCM(){}

public:
	BYTE tag[16];  // validation tag for gcm
};

// Todo: Keychain class: object[blockID][baseObjects] or object[blockID].atINFO(id) or object[blockID].atFORMAT(id)
// Keychain write: blockID, n_info, n_format | baseID : key, iv | ...
struct Keychainchain{
public:
	typedef Keychainchain     self_type;
	typedef KeychainEntryGCM  value_type;
	typedef value_type&       reference;
	typedef const value_type& const_reference;

public:
	Keychainchain() : blockID(0), info_keys(nullptr), format_keys(nullptr){}
	Keychainchain(const U32 n_info_fields, const U32 n_format_fields) : blockID(0), info_keys(new value_type[n_info_fields]), format_keys(new value_type[n_format_fields]){}
	~Keychainchain(){
		delete [] this->info_keys;
		delete [] this->format_keys;
	}
	// Todo: rule of 5

	inline U64& getBlockID(void){ return(this->blockID); }
	inline const U64& getBlockID(void) const{ return(this->blockID); }

	inline reference operator[](const U32 position){ return(this->base_keys[position]); }
	inline const_reference operator[](const U32 position) const{ return(this->base_keys[position]); }
	inline reference at(const U32 position){ return(this->base_keys[position]); }
	inline const_reference at(const U32 position) const{ return(this->base_keys[position]); }
	inline reference atINFO(const U32 position){ return(this->info_keys[position]); }
	inline const_reference atINFO(const U32 position) const{ return(this->info_keys[position]); }
	inline reference atFORMAT(const U32 position){ return(this->format_keys[position]); }
	inline const_reference atFORMAT(const U32 position) const{ return(this->format_keys[position]); }

public:
	U64 blockID;
	value_type  base_keys[19];
	value_type* info_keys;
	value_type* format_keys;
};

class Keychain{
private:
    typedef Keychain           self_type;
	typedef std::size_t        size_type;
    typedef Keychainchain      value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

public:
    Keychain();
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

private:
    size_type n_entries_;
    size_type n_capacity_;
    pointer   entries_;
	// todo: hash table for blockID
};

class EncryptionDecorator{
public:
	typedef EncryptionDecorator       self_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer stream_container;
	typedef io::BasicBuffer           buffer_type;
	typedef KeychainEntryGCM          value_type;

public:
	bool encryptAES256(variant_block_type& block){
		Keychainchain keychain(block.footer.n_info_streams, block.footer.n_format_streams);

		//value_type entry;
		if(!this->encryptAES256(block.meta_contig_container,keychain[0])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_contig_container, keychain[0]);

		if(!this->encryptAES256(block.meta_positions_container,keychain[1])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_positions_container, keychain[1]);

		if(!this->encryptAES256(block.meta_refalt_container,keychain[2])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_refalt_container, keychain[2]);

		if(!this->encryptAES256(block.meta_controller_container,keychain[3])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_controller_container, keychain[3]);

		if(!this->encryptAES256(block.meta_quality_container,keychain[4])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_quality_container, keychain[4]);

		if(!this->encryptAES256(block.meta_names_container,keychain[5])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_names_container, keychain[5]);

		if(!this->encryptAES256(block.gt_rle8_container,keychain[6])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_rle8_container, keychain[6]);

		if(!this->encryptAES256(block.gt_rle16_container,keychain[7])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_rle16_container, keychain[7]);

		if(!this->encryptAES256(block.gt_rle32_container,keychain[8])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_rle32_container, keychain[8]);

		if(!this->encryptAES256(block.gt_rle64_container,keychain[9])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_rle64_container, keychain[9]);

		if(!this->encryptAES256(block.meta_alleles_container,keychain[10])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_alleles_container, keychain[10]);

		if(!this->encryptAES256(block.gt_simple8_container,keychain[11])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_simple8_container, keychain[11]);

		if(!this->encryptAES256(block.gt_simple16_container,keychain[12])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_simple16_container, keychain[12]);

		if(!this->encryptAES256(block.gt_simple32_container,keychain[13])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_simple32_container, keychain[13]);

		if(!this->encryptAES256(block.gt_simple64_container,keychain[14])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_simple64_container, keychain[14]);

		if(!this->encryptAES256(block.gt_support_data_container,keychain[15])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.gt_support_data_container, keychain[15]);

		if(!this->encryptAES256(block.meta_info_map_ids,keychain[16])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_info_map_ids, keychain[16]);

		if(!this->encryptAES256(block.meta_filter_map_ids,keychain[17])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_filter_map_ids, keychain[17]);

		if(!this->encryptAES256(block.meta_format_map_ids,keychain[18])){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
		this->decryptAES256(block.meta_format_map_ids, keychain[18]);

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(!this->encryptAES256(block.info_containers[i], keychain.atINFO(i))){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
			this->decryptAES256(block.info_containers[i], keychain.atINFO(i));
		}
		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
			if(!this->encryptAES256(block.format_containers[i], keychain.atFORMAT(i))){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
			this->decryptAES256(block.format_containers[i], keychain.atFORMAT(i));
		}

		return(true);
	}

	bool encryptAES256(stream_container& container, value_type& entry){
		// temp
		entry.initiateRandom();

		EVP_CIPHER_CTX *ctx = NULL;
		int len = 0;

		/* Create and initialise the context */
		if(!(ctx = EVP_CIPHER_CTX_new())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
			return false;
		}

		/* Initialise the encryption operation. */
		if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption..." << std::endl;
			return false;
		}

		// 16 bytes = 128 bits
		if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption tag..." << std::endl;
			return false;
		}

		/* Initialise key and IV */
		if(1 != EVP_EncryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption key and IV..." << std::endl;
			return false;
		}

		/* Provide the message to be encrypted, and obtain the encrypted output.
		 * EVP_EncryptUpdate can be called multiple times if necessary
		 */
		this->buffer.reset();
		this->buffer.resize(container.buffer_data.size() + container.buffer_strides.size() + 65536);
		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)this->buffer.data(), &len, (BYTE*)container.buffer_data.data(), container.buffer_data.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}
		this->buffer.n_chars = len;

		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)&this->buffer[this->buffer.n_chars], &len, (BYTE*)container.buffer_strides.data(), container.buffer_strides.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		/* Finalise the encryption. Normally ciphertext bytes may be written at
		 * this stage, but this does not occur in GCM mode
		 */
		if(1 != EVP_EncryptFinal_ex(ctx, (BYTE*)this->buffer.data() + len, &len)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to finalise the encryption..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		/* Get the tag */
		container.buffer_data.resize(this->buffer.size());
		if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, &entry.tag[0])){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to retrieve the GCM tag..." << std::endl;
			return false;
		}

		/* Clean up */
		EVP_CIPHER_CTX_free(ctx);

		// Trigger encryption flag
		container.header.data_header.controller.encryption = YON_ENCRYPTION_AES_256_GCM;
		memcpy(container.buffer_data.data(), this->buffer.data(), this->buffer.size());
		container.buffer_data.n_chars = this->buffer.size();
		container.header.data_header.eLength = this->buffer.size();
		//std::cerr << container.header.data_header.eLength << "/" << this->buffer.size() << '\t' << container.header.data_header.cLength << "/" << container.header.stride_header.cLength << std::endl;
		assert(container.header.data_header.eLength == container.header.data_header.cLength + container.header.stride_header.cLength);

		return(true);
	}

	bool decryptAES256(stream_container& container, value_type& entry){
		if(container.header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
			return true;

		if(container.header.data_header.controller.encryption != YON_ENCRYPTION_AES_256_GCM){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Wrong decryption algorithm used..." << std::endl;
			return false;
		}

		EVP_CIPHER_CTX *ctx = NULL;
		int len = 0, plaintext_len = 0, ret;

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

		/* Set IV length. Not necessary if this is 12 bytes (96 bits) */
		if(!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the IV context..." << std::endl;
			return(false);
		}

		/* Initialise key and IV */
		if(!EVP_DecryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the decryption..." << std::endl;
			return(false);
		}

		/* Provide the message to be decrypted, and obtain the plaintext output.
		 * EVP_DecryptUpdate can be called multiple times if necessary
		 */
		this->buffer.reset();
		if(container.buffer_data.size()){
			//std::cerr << "decrypting: " << container.buffer_data.size() << std::endl;
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
			/* Success */
			plaintext_len += len;
			container.buffer_data.reset();
			container.buffer_strides.reset();
			container.buffer_data.resize(container.header.data_header.cLength + 65536);
			container.buffer_strides.resize(container.header.stride_header.cLength + 65536);
			memcpy(container.buffer_data.data(), this->buffer.data(), container.header.data_header.cLength);
			memcpy(container.buffer_strides.data(), &this->buffer[container.header.data_header.cLength], container.header.stride_header.cLength);
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
