#ifndef CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_
#define CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_

#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <cassert>

#include "io/basic_buffer.h"

namespace tachyon{
namespace containers{

// Controller type for stream container
struct DataContainerHeaderController{
public:
	typedef DataContainerHeaderController self_type;
	typedef io::BasicBuffer               buffer_type;

public:
	DataContainerHeaderController() :
		signedness(0),
		mixedStride(0),
		type(0),
		encoder(0),
		uniform(0),
		encryption(0),
		preprocessor(0)
	{}

	~DataContainerHeaderController(){}

	inline void clear(){
		this->signedness   = 0;
		this->mixedStride  = 0;
		this->type         = 0;
		this->encoder      = 0;
		this->uniform      = 0;
		this->encryption   = 0;
		this->preprocessor = 0;
	}

	inline bool IsEncrypted(void) const{ return(this->encryption > 0); }
	inline bool CompareType(const uint8_t& type) const{ return(this->type == type); }
	inline bool CompareTypeSign(const uint8_t& type, const bool& sign) const{ return(this->type == type && this->signedness == sign); }

	self_type& operator=(const self_type& other){
		this->signedness   = other.signedness;
		this->mixedStride  = other.mixedStride;
		this->type         = other.type;
		this->encoder      = other.encoder;
		this->uniform      = other.uniform;
		this->encryption   = other.encryption;
		this->preprocessor = other.preprocessor;
		return(*this);
	}

	bool operator==(const self_type& other) const{
		if(this->signedness   != other.signedness)   return false;
		if(this->mixedStride  != other.mixedStride)  return false;
		if(this->type         != other.type)         return false;
		if(this->encoder      != other.encoder)      return false;
		if(this->uniform      != other.uniform)      return false;
		if(this->encryption   != other.encryption)   return false;
		if(this->preprocessor != other.preprocessor) return false;
		return true;
	}
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

private:
	friend buffer_type& operator<<(buffer_type& buffer,const self_type& controller){
		const uint32_t c = controller.signedness << 0  |
					  controller.mixedStride  << 1  |
					  controller.type         << 2  |
					  controller.encoder      << 8  |
					  controller.uniform      << 13 |
					  controller.encryption   << 14 |
					  controller.preprocessor << 16;

		//const uint16_t* c = reinterpret_cast<const uint16_t* const>(&controller);
		buffer += c;
		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller){
		const uint32_t c = controller.signedness  << 0  |
					  controller.mixedStride  << 1  |
					  controller.type         << 2  |
					  controller.encoder      << 8  |
					  controller.uniform      << 13 |
					  controller.encryption   << 14 |
					  controller.preprocessor << 16;

		//assert(*reinterpret_cast<const uint16_t* const>(&controller) == c);

		stream.write(reinterpret_cast<const char*>(&c), sizeof(uint32_t));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		stream.read(reinterpret_cast<char*>(&controller), sizeof(uint32_t));
		return(stream);
	}

	friend buffer_type& operator>>(buffer_type& buffer, self_type& controller){
		uint32_t* c = reinterpret_cast<uint32_t*>(&controller);
		buffer >> *c;
		return(buffer);
	}

public:
	uint32_t signedness:  1, // Signed type
	         mixedStride: 1, // Different stride sizes
	         type:        6, // Base typing (extra bits reserved for future use)
	         encoder:     5, // Encoder bits (see encoder for values)
	         uniform:     1, // Triggered if all values in the buffer are the same
	         encryption:  2, // Encryption type
			 preprocessor: 16; // preprocessor bits (extra reserved for future used)
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
