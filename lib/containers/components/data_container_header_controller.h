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

	inline bool isEncrypted(void) const{ return(this->encryption > 0); }
	inline bool compareType(const uint8_t& type) const{ return(this->type == type); }
	inline bool compareTypeSign(const uint8_t& type, const bool& sign) const{ return(this->type == type && this->signedness == sign); }

	self_type& operator=(const self_type& other){
		this->signedness  = other.signedness;
		this->mixedStride = other.mixedStride;
		this->type        = other.type;
		this->encoder     = other.encoder;
		this->uniform     = other.uniform;
		this->encryption  = other.encryption;
		return(*this);
	}

	bool operator==(const self_type& other) const{
		if(this->signedness  != other.signedness)  return false;
		if(this->mixedStride != other.mixedStride) return false;
		if(this->type        != other.type)        return false;
		if(this->encoder     != other.encoder)     return false;
		if(this->uniform     != other.uniform)     return false;
		if(this->encryption  != other.encryption)  return false;
		return true;
	}
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

private:
	friend buffer_type& operator<<(buffer_type& buffer,const self_type& controller){
		const uint16_t c = controller.signedness  << 0  |
					  controller.mixedStride << 1  |
					  controller.type        << 2  |
					  controller.encoder     << 8  |
					  controller.uniform     << 13 |
					  controller.encryption  << 14;

		//const uint16_t* c = reinterpret_cast<const uint16_t* const>(&controller);
		buffer += c;
		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller){
		const uint16_t c = controller.signedness  << 0  |
					  controller.mixedStride << 1  |
					  controller.type        << 2  |
					  controller.encoder     << 8  |
					  controller.uniform     << 13 |
					  controller.encryption  << 14;

		//assert(*reinterpret_cast<const uint16_t* const>(&controller) == c);

		stream.write(reinterpret_cast<const char*>(&c), sizeof(uint16_t));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		stream.read(reinterpret_cast<char*>(&controller), sizeof(uint16_t));
		return(stream);
	}

	friend buffer_type& operator>>(buffer_type& buffer, self_type& controller){
		uint16_t* c = reinterpret_cast<uint16_t*>(&controller);
		buffer >> *c;
		return(buffer);
	}

public:
	uint16_t signedness:  1, // Signed type
		mixedStride: 1, // Different stride sizes
		type:        6, // Base typing (extra bits saved for future use)
		encoder:     5, // Encoder bits (see encoder for values)
		uniform:     1, // Triggered if all values in the buffer are the same
		encryption:  2; // Encryption type
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
