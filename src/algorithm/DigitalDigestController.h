#ifndef ALGORITHM_DIGITALDIGESTCONTROLLER_H_
#define ALGORITHM_DIGITALDIGESTCONTROLLER_H_

#include <openssl/sha.h>

namespace Tachyon{
namespace Algorithm{

class DigitalDigestController{
private:
	typedef DigitalDigestController self_type;

public:
	DigitalDigestController(void){}
	~DigitalDigestController(void){}

	inline bool initialize(){
		if(!SHA512_Init(&this->context))
			return false;

		return true;
	}

	inline bool update(void* input, const U32 length){
		if(!SHA512_Update(&this->context, (BYTE*)input, length))
			return false;

		return true;
	}

	inline bool finalize(){
		if(!SHA512_Final(&this->sha512_digest[0], &context))
			return false;

		return true;
	}

	inline void clear(void){ memset(&this->sha512_digest[0], 0, 64); }

public:
	SHA512_CTX context;
	BYTE sha512_digest[64];
};

}
}

#endif /* ALGORITHM_DIGITALDIGESTCONTROLLER_H_ */
