#ifndef ALGORITHM_ENCRYPTION_ENCRYPTIONENTRY_H_
#define ALGORITHM_ENCRYPTION_ENCRYPTIONENTRY_H_

namespace Tachyon{
namespace Encryption{

struct EncryptionEntry{
	typedef EncryptionEntry self_type;

public:
	BYTE controller;
	BYTE iv[16]; // 128 bit initiation vector
};

}
}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTIONENTRY_H_ */
