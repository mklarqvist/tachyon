#ifndef CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_
#define CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include "../../support/type_definitions.h"

namespace tachyon{
namespace containers{

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

	inline void clear(){
		this->signedness  = 0;
		this->mixedStride = 0;
		this->type        = 0;
		this->encoder     = 0;
		this->uniform     = 0;
		this->encryption  = 0;
	}

	inline const bool isEncrypted(void) const{ return(this->encryption > 0); }
	inline const bool compareType(const BYTE& type) const{ return(this->type == type); }
	inline const bool compareTypeSign(const BYTE& type, const bool& sign) const{ return(this->type == type && this->signedness == sign); }

	self_type& operator=(const self_type& other){
		this->signedness  = other.signedness;
		this->mixedStride = other.mixedStride;
		this->type        = other.type;
		this->encoder     = other.encoder;
		this->uniform     = other.uniform;
		this->encryption  = other.encryption;
		return(*this);
	}

	const bool operator==(const self_type& other) const{
		if(this->signedness  != other.signedness)  return false;
		if(this->mixedStride != other.mixedStride) return false;
		if(this->type        != other.type)        return false;
		if(this->encoder     != other.encoder)     return false;
		if(this->uniform     != other.uniform)     return false;
		if(this->encryption  != other.encryption)  return false;
		return true;
	}
	inline const bool operator!=(const self_type& other) const{ return(!(*this == other)); }

private:
	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer,const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		buffer += *c;
		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller){
		const U16 c = controller.signedness  << 0  |
					  controller.mixedStride << 1  |
					  controller.type        << 2  |
					  controller.encoder     << 8  |
					  controller.uniform     << 13 |
					  controller.encryption  << 14;

		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& controller){
		stream.read(reinterpret_cast<char*>(&controller), sizeof(U16));
		return(stream);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		*c << buffer;
		return(buffer);
	}

public:
	U16 signedness:  1, // Signed type
		mixedStride: 1, // Different stride sizes
		type:        6, // Base typing (extra bits saved for future use)
		encoder:     5, // Encoder bits (see encoder for values)
		uniform:     1, // Triggered if all values in the buffer are the same
		encryption:  2; // Encryption type
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
