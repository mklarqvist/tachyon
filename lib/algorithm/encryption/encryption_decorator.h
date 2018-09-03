#ifndef ALGORITHM_ENCRYPTION_ENCRYPTION_DECORATOR_H_
#define ALGORITHM_ENCRYPTION_ENCRYPTION_DECORATOR_H_

#include <openssl/evp.h>
#include <openssl/aes.h>
#include <openssl/rand.h>

#include "keychain.h"
#include "containers/variant_block.h"
#include "support/enums.h"


namespace tachyon{

class EncryptionDecorator{
public:
	typedef EncryptionDecorator       self_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer stream_container;
	typedef io::BasicBuffer           buffer_type;
	typedef Keychain                  keychain_type;

public:
	EncryptionDecorator() = default;
	~EncryptionDecorator() = default;
	bool Encrypt(variant_block_type& block, keychain_type& keychain, TACHYON_ENCRYPTION encryption_type);
	bool Decrypt(variant_block_type& block, keychain_type& keychain);
	bool EncryptAES256(variant_block_type& block, keychain_type& keychain);
	bool EncryptAES256(stream_container& container, keychain_type& keychain);
	bool DecryptAES256(stream_container& container, keychain_type& keychain);

public:
	buffer_type buffer;
};

}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTION_DECORATOR_H_ */
