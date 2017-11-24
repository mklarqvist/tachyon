#ifndef ALGORITHM_ENCRYPTION_OPENSSL_AES_H_
#define ALGORITHM_ENCRYPTION_OPENSSL_AES_H_

#include <openssl/evp.h>
#include <openssl/rand.h>

namespace Tachyon{
namespace Encryption{

int aes_init(BYTE* key_data, int key_data_len, BYTE* retKey, BYTE* retIV){
  int i, nrounds = 14;

  /*
   * Gen key & IV for AES 256 CBC mode. A SHA1 digest is used to hash the supplied key material.
   * nrounds is the number of times the we hash the material. More rounds are more secure but
   * slower.
   */
  i = EVP_BytesToKey(EVP_aes_256_cbc(), EVP_sha256(), nullptr, key_data, key_data_len, nrounds, retKey, retIV);
  if (i != 32) {
    printf("Key size is %d bits - should be 256 bits\n", i);
    return -1;
  }

  return 0;
}

inline int aes_encrypt(BYTE* input, S32 length, BYTE* key, BYTE* iv, BYTE* output){
	EVP_CIPHER_CTX *ctx;
	int len;
	int ciphertext_len;

	/* Create and initialise the context */
	if(!(ctx = EVP_CIPHER_CTX_new())){
		exit(1);
	}

	/* Initialise the encryption operation. IMPORTANT - ensure you use a key
	* and IV size appropriate for your cipher
	* In this example we are using 256 bit AES (i.e. a 256 bit key). The
	* IV size for *most* modes is the same as the block size. For AES this
	* is 128 bits */
	if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv)){
		exit(1);
	}

	/* Provide the message to be encrypted, and obtain the encrypted output.
	* EVP_EncryptUpdate can be called multiple times if necessary
	*/
	if(1 != EVP_EncryptUpdate(ctx, output, &len, input, length)){
		exit(1);
	}
	ciphertext_len = len;

	/* Finalise the encryption. Further ciphertext bytes may be written at
	* this stage.
	*/
	if(1 != EVP_EncryptFinal_ex(ctx, output + len, &len)){
		exit(1);
	}
	ciphertext_len += len;

	/* Clean up */
	EVP_CIPHER_CTX_free(ctx);

	return ciphertext_len;
}

inline int aes_decrypt(BYTE* input, S32 length, BYTE* key, BYTE* iv, BYTE* output){
	EVP_CIPHER_CTX *ctx;

	int len;

	int plaintext_len;

	/* Create and initialise the context */
	if(!(ctx = EVP_CIPHER_CTX_new())){
		exit(1);
	}

	/* Initialise the decryption operation. IMPORTANT - ensure you use a key
	* and IV size appropriate for your cipher
	* In this example we are using 256 bit AES (i.e. a 256 bit key). The
	* IV size for *most* modes is the same as the block size. For AES this
	* is 128 bits */
	if(1 != EVP_DecryptInit_ex(ctx, EVP_aes_256_cbc(), NULL, key, iv)){
		exit(1);
	}

	/* Provide the message to be decrypted, and obtain the plaintext output.
	* EVP_DecryptUpdate can be called multiple times if necessary
	*/
	if(1 != EVP_DecryptUpdate(ctx, output, &len, input, length)){
		std::cerr << "update fail" << std::endl;
		exit(1);
	}
	plaintext_len = len;

	/* Finalise the decryption. Further plaintext bytes may be written at
	* this stage.
	*/
	if(1 != EVP_DecryptFinal_ex(ctx, output + len, &len)){
		std::cerr << "final fail" << std::endl;
		//exit(1);
	}
	plaintext_len += len;

	/* Clean up */
	EVP_CIPHER_CTX_free(ctx);

	return plaintext_len;
}

}
}
#endif /* ALGORITHM_ENCRYPTION_OPENSSL_AES_H_ */
