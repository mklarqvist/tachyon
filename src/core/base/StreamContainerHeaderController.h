#ifndef CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_
#define CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include "../../support/TypeDefinitions.h"

namespace Tachyon{
namespace Core{

enum CORE_ENCRYPTION{ENCRYPTION_NONE, ENCRYPTION_AES_128, ENCRYPTION_AES_256, ENCRYPTION_RSA_4096};
enum CORE_ENCODING{ENCODING_NONE, ENCODING_DEFLATE};

// Controller type for stream container
struct StreamContainerHeaderController{
	typedef StreamContainerHeaderController self_type;

public:
	StreamContainerHeaderController() :
		signedness(0),
		mixedStride(0),
		type(0),
		encoder(0),
		uniform(0),
		encryption(0)
	{}
	~StreamContainerHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }
	inline bool isEncrypted(void) const{ return(this->encryption > 0); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		stream.write(reinterpret_cast<const char*>(c), sizeof(U16));
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& controller){
		stream.read(reinterpret_cast<char*>(&controller), sizeof(U16));
		return(stream);
	}

public:
	// 6 base values (4 integers + 2 floats)
	U16 signedness: 1,
		mixedStride: 1,
		type: 6,    // base typing (extra bits saved for future use)
		encoder: 5, // encoder bits (0 = uncompressed)
		uniform: 1, // triggered if all values are the same
		encryption: 2;
};

}
}



#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
