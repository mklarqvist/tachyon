#ifndef ALGORITHM_ENCRYPTION_OPENSSL_AES_H_
#define ALGORITHM_ENCRYPTION_OPENSSL_AES_H_

#include <openssl/evp.h>

int aes_init(unsigned char *key_data, int key_data_len, unsigned char *salt, EVP_CIPHER_CTX *e_ctx, EVP_CIPHER_CTX *d_ctx);
unsigned char *aes_encrypt(EVP_CIPHER_CTX *e, unsigned char *plaintext, int *len);
unsigned char *aes_decrypt(EVP_CIPHER_CTX *e, unsigned char *ciphertext, int *len);



#endif /* ALGORITHM_ENCRYPTION_OPENSSL_AES_H_ */
