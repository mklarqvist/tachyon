#ifndef CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_
#define CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include "../../support/type_definitions.h"

namespace tachyon{
namespace containers{

enum TACHYON_ENCRYPTION{YON_ENCRYPTION_NONE,
	                    YON_ENCRYPTION_AES_128,
					    YON_ENCRYPTION_AES_256,
					    YON_ENCRYPTION_RSA_4096};

// Controller type for stream container
struct DataContainerHeaderController{
	typedef DataContainerHeaderController self_type;

public:
	DataContainerHeaderController() :
		signedness(0),
		mixedStride(0),
		type(0),
		encoder(0),
		uniform(0),
		encryption(0)
	{}
	~DataContainerHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }
	inline bool isEncrypted(void) const{ return(this->encryption > 0); }

	inline const bool compareType(const BYTE& type) const{ return(this->type == type); }
	inline const bool compareTypeSign(const BYTE& type, const bool& sign) const{ return(this->type == type && this->signedness == sign); }

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer,const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		buffer += *c;
		return(buffer);
	}

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
	U16 signedness: 1, // Signed type
		mixedStride: 1,// Different stride sizes
		type: 6,       // Base typing (extra bits saved for future use)
		encoder: 5,    // Encoder bits (see encoder for values)
		uniform: 1,    // Triggered if all values in the buffer are the same
		encryption: 2; // Encryption type
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
