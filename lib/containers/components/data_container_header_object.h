#ifndef CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_
#define CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_

#include "io/basic_buffer.h"
#include "support/enums.h"
#include "data_container_header_controller.h"

namespace tachyon{
namespace containers{

struct DataContainerHeaderObject{
	typedef DataContainerHeaderObject     self_type;
	typedef DataContainerHeaderController controller_type;

	DataContainerHeaderObject();
	DataContainerHeaderObject(const self_type& other);
	DataContainerHeaderObject(self_type&& other) noexcept;
	self_type& operator=(const self_type& other);
	self_type& operator=(self_type&& other) noexcept;
	~DataContainerHeaderObject();

	void reset(void);
	const bool operator==(const self_type& other) const;
	inline const bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	const SBYTE getPrimitiveWidth(void) const;

	//
	inline S32& getStride(void){ return(this->stride); }
	inline const S32& getStride(void) const{ return(this->stride); }

	inline const bool isUniform(void) const{ return(this->controller.uniform); }
	inline const bool isSigned(void) const{ return(this->controller.signedness); }
	inline const bool hasMixedStride(void) const{ return(this->controller.mixedStride); }
	inline void setUniform(const bool yes){ this->controller.uniform = yes; }
	inline void setSignedness(const bool yes){ this->controller.signedness = yes; }
	inline void setMixedStride(const bool yes){ this->controller.mixedStride = yes; }

	inline const TACHYON_CORE_TYPE getPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->controller.type)); }
	inline const TACHYON_CORE_COMPRESSION getEncoder(void) const{ return(TACHYON_CORE_COMPRESSION(this->controller.encoder)); }

	// Set types
	inline void setType(const TACHYON_CORE_TYPE& type){ this->controller.type = type; }

	// Checksum
	inline U32& getChecksum(void){ return(this->crc); }
	inline const U32& getChecksum(void) const{ return(this->crc); }
	inline const bool checkChecksum(const U32 checksum) const{ return(this->crc == checksum); }

private:
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry){
		buffer << entry.controller;
		buffer += entry.stride;
		buffer += entry.offset;
		buffer += entry.cLength;
		buffer += entry.uLength;
		buffer += entry.eLength;
		buffer += entry.crc;
		buffer += entry.global_key;
		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.stride),    sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&entry.offset),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.cLength),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uLength),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.eLength),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.crc),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.global_key),sizeof(S32));
		return(stream);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& entry){
		buffer >> entry.controller;
		buffer >> entry.stride;
		buffer >> entry.offset;
		buffer >> entry.cLength;
		buffer >> entry.uLength;
		buffer >> entry.eLength;
		buffer >> entry.crc;
		buffer >> entry.global_key;
		return(buffer);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.stride),     sizeof(S32));
		stream.read(reinterpret_cast<char*>(&entry.offset),     sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.cLength),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uLength),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.eLength),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.crc),        sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.global_key), sizeof(S32));

		return(stream);
	}

public:
	controller_type controller; // controller bits
	S32 stride;                 // stride size: -1 if not uniform, a non-zero positive value otherwise
	U32 offset;                 // relative file offset
	U32 cLength;                // compressed length
	U32 uLength;                // uncompressed length
	U32 eLength;                // encrypted length
	U32 crc;                    // crc32 checksum
	S32 global_key;             // global key
};

}
}



#endif /* CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_ */
