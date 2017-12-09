#ifndef ALGORITHM_DIGITALDIGESTCONTROLLER_H_
#define ALGORITHM_DIGITALDIGESTCONTROLLER_H_

#include <openssl/sha.h>

namespace Tachyon{
namespace Algorithm{

class DigitalDigestController{
private:
	typedef DigitalDigestController self_type;
	typedef Core::StreamContainer container_type;

public:
	DigitalDigestController(void){}
	~DigitalDigestController(void){}

	inline bool initialize(){
		if(!SHA512_Init(&this->context))
			return false;

		return true;
	}

	inline bool update(const container_type& container){
		if(!this->update(container.buffer_data.data, container.buffer_data.pointer)){
			std::cerr << "failed update" << std::endl;
			return false;
		}
		if(container.header.controller.mixedStride){
			if(!this->update(container.buffer_strides.data, container.buffer_strides.pointer)){
				std::cerr << "failed update" << std::endl;
				return false;
			}
		}
		return true;
	}

	inline bool updateUncompressed(const container_type& container){
		if(!this->update(container.buffer_data_uncompressed.data, container.buffer_data_uncompressed.pointer)){
			std::cerr << "failed update" << std::endl;
			return false;
		}
		if(container.header.controller.mixedStride){
			if(!this->update(container.buffer_strides_uncompressed.data, container.buffer_strides_uncompressed.pointer)){
				std::cerr << "failed update" << std::endl;
				return false;
			}
		}
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

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.sha512_digest), 64);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.sha512_digest), 64);
		return(stream);
	}

public:
	SHA512_CTX context;
	BYTE sha512_digest[64];
};

}
}

#endif /* ALGORITHM_DIGITALDIGESTCONTROLLER_H_ */
