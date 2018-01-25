#ifndef ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_
#define ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_

namespace tachyon{
namespace encryption{

struct EncryptionDecorator{
	typedef EncryptionDecorator self_type;
	typedef core::DataContainer stream_container;
	typedef io::BasicBuffer buffer_type;

	bool encrypt(stream_container& container);
	bool decrypt(stream_container& container);

public:
	BYTE key[32];  // 256 bit key
	BYTE iv[16];   // 128 bit initiation vector
	BYTE salt[12]; // salt for taste
	buffer_type buffer;
};

}
}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTIONDECORATOR_H_ */
