#include "encryption_decorator.h"

namespace tachyon{
namespace encryption{

bool EncryptionDecorator::encrypt(variant_block_type& block, keychain_type& keychain, TACHYON_ENCRYPTION encryption_type){
	if(encryption_type == YON_ENCRYPTION_AES_256_GCM){
		return(this->encryptAES256(block, keychain));
	}
	return(false);
}

bool EncryptionDecorator::decryptAES256(variant_block_type& block, keychain_type& keychain){

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(!this->decryptAES256(block.base_containers[i], keychain)){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl;
			return false;
		}
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(!this->decryptAES256(block.info_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(!this->decryptAES256(block.format_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to decrypt!" << std::endl; return false; }
	}

	return(true);
}

bool EncryptionDecorator::encryptAES256(variant_block_type& block, keychain_type& keychain){
	uint8_t RANDOM_BYTES[32];
	RAND_bytes(&RANDOM_BYTES[0], 32);
	block.header.block_hash = XXH64(&RANDOM_BYTES[0], 32, 1337);

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(!this->encryptAES256(block.base_containers[i], keychain)){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl;
			return false;
		}
	}
	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(!this->encryptAES256(block.info_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(!this->encryptAES256(block.format_containers[i], keychain)){ std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl; return false; }
	}

	return(true);
}

bool EncryptionDecorator::encryptAES256(stream_container& container, keychain_type& keychain){
	aes256gcm_type entry;
	entry.initiateRandom();
	entry.encryption_type = YON_ENCRYPTION_AES_256_GCM;

	EVP_CIPHER_CTX *ctx = NULL;
	int32_t len = 0;

	if(!(ctx = EVP_CIPHER_CTX_new())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
		return false;
	}

	if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption..." << std::endl;
		return false;
	}

	// 16 uint8_ts = 128 bits
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

	if(1 != EVP_EncryptUpdate(ctx, (uint8_t*)this->buffer.data(), &len, (uint8_t*)temp.data(), temp.size())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars = len;

	if(1 != EVP_EncryptUpdate(ctx, (uint8_t*)&this->buffer[this->buffer.n_chars], &len, (uint8_t*)container.buffer_data.data(), container.buffer_data.size())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars += len;

	if(1 != EVP_EncryptUpdate(ctx, (uint8_t*)&this->buffer[this->buffer.n_chars], &len, (uint8_t*)container.buffer_strides.data(), container.buffer_strides.size())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars += len;

	if(1 != EVP_EncryptFinal_ex(ctx, (uint8_t*)this->buffer.data() + len, &len)){
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

	const uint64_t hashID = keychain.getRandomHashIdentifier();
	container.header.identifier = hashID;
	entry.field_id = hashID;
	keychain += entry; // add key to keychain

	return(true);
}

bool EncryptionDecorator::decryptAES256(stream_container& container, keychain_type& keychain){
	if(container.buffer_data.size() == 0)
		return true;

	if(container.header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
		return true;

	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_AES_256_GCM){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Wrong decryption algorithm used..." << std::endl;
		return false;
	}

	uint32_t* match = nullptr;
	if(keychain.getHashIdentifier(container.header.identifier, match) == false){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Did not find ID in keychain..." << std::endl;
		return false;
	}
	KeychainKeyGCM<>& entry = keychain[*match];

	EVP_CIPHER_CTX *ctx = NULL;
	int32_t len = 0, plaintext_len = 0, ret;

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
		if(!EVP_DecryptUpdate(ctx, (uint8_t*)this->buffer.data(), &len, (uint8_t*)container.buffer_data.data(), container.buffer_data.size())){
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
	ret = EVP_DecryptFinal_ex(ctx, (uint8_t*)this->buffer.data() + len, &len);
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

}
}
