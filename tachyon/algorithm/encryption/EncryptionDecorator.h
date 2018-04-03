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
	~KeychainEntry(){}

	void print(void){
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

// Todo: Keychain class: object[blockID][baseObjects] or object[blockID].atINFO(id) or object[blockID].atFORMAT(id)
// Keychain write: blockID, n_info, n_format | baseID : key, iv | ...

struct EncryptionDecorator{
public:
	typedef EncryptionDecorator       self_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer stream_container;
	typedef io::BasicBuffer           buffer_type;
	typedef KeychainEntry             value_type;

public:
	bool encrypt(variant_block_type& block){
		this->encrypt(block.meta_contig_container);
		this->encrypt(block.meta_positions_container);
		this->encrypt(block.meta_refalt_container);
		this->encrypt(block.meta_controller_container);
		this->encrypt(block.meta_quality_container);
		this->encrypt(block.meta_names_container);
		this->encrypt(block.gt_rle8_container);
		this->encrypt(block.gt_rle16_container);
		this->encrypt(block.gt_rle32_container);
		this->encrypt(block.gt_rle64_container);
		this->encrypt(block.meta_alleles_container);
		this->encrypt(block.gt_simple8_container);
		this->encrypt(block.gt_simple16_container);
		this->encrypt(block.gt_simple32_container);
		this->encrypt(block.gt_simple64_container);
		this->encrypt(block.gt_support_data_container);
		this->encrypt(block.meta_info_map_ids);
		this->encrypt(block.meta_filter_map_ids);
		this->encrypt(block.meta_format_map_ids);
		for(U32 i = 0; i < block.footer.n_info_streams; ++i)   this->encrypt(block.info_containers[i]);
		for(U32 i = 0; i < block.footer.n_format_streams; ++i) this->encrypt(block.format_containers[i]);

		return(true);
	}

	bool encrypt(stream_container& container){
		// Todo:
		// 1) Generate N random bytes
		// 2) encrypt N + N_data + data
		// 3) store
		value_type entry;
		entry.initiateRandom();
		//entry.print();
		//std::cout.write((char*)&entry.key[0], 16);
		//std::cout.write((char*)&entry.iv[0],  12);


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
		this->buffer.resize(container.buffer_data.size() + 65536);
		if(1 != EVP_EncryptUpdate(ctx, (BYTE*)this->buffer.data(), &len, (BYTE*)container.buffer_data.data(), container.buffer_data.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
			return false;
		}

		this->buffer.n_chars = len;

		/* Finalise the encryption. Normally ciphertext bytes may be written at
		 * this stage, but this does not occur in GCM mode
		 */
		if(1 != EVP_EncryptFinal_ex(ctx, (BYTE*)this->buffer.data() + len, &len)){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to finalise the encryption..." << std::endl;
			return false;
		}
		this->buffer.n_chars += len;

		/* Get the tag */
		//if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag))
		//	exit(1);

		/* Clean up */
		EVP_CIPHER_CTX_free(ctx);

		// Trigger encryption flag
		container.header.data_header.controller.encryption = YON_ENCRYPTION_AES_256;
		memcpy(container.buffer_data.data(), this->buffer.data(), this->buffer.size());
		assert(container.header.data_header.cLength == this->buffer.size());
		return(true);
	}

	bool decrypt(stream_container& container);

public:
	buffer_type buffer;
};

}
}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_ */
