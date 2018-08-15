#ifndef CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_
#define CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_

#include "io/basic_buffer.h"
#include "support/enums.h"
#include "data_container_header_controller.h"

#include "openssl/md5.h"

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
	bool operator==(const self_type& other) const;
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	SBYTE GetPrimitiveWidth(void) const;

	//
	inline S32& GetStride(void){ return(this->stride); }
	inline const S32& GetStride(void) const{ return(this->stride); }

	inline bool IsUniform(void) const{ return(this->controller.uniform); }
	inline bool IsSigned(void) const{ return(this->controller.signedness); }
	inline bool HasMixedStride(void) const{ return(this->controller.mixedStride); }
	inline void SetUniform(const bool yes){ this->controller.uniform = yes; }
	inline void SetSignedness(const bool yes){ this->controller.signedness = yes; }
	inline void SetMixedStride(const bool yes){ this->controller.mixedStride = yes; }

	inline TACHYON_CORE_TYPE GetPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->controller.type)); }
	inline TACHYON_CORE_COMPRESSION GetEncoder(void) const{ return(TACHYON_CORE_COMPRESSION(this->controller.encoder)); }

	// Set types
	inline void SetType(const TACHYON_CORE_TYPE& type){ this->controller.type = type; }

	// Checksum
	inline uint8_t* GetChecksum(void){ return(&this->crc[0]); }
	inline const uint8_t* GetChecksum(void) const{ return(&this->crc[0]); }
	bool CheckChecksum(const uint8_t* compare) const{
		for(U32 i = 0; i < MD5_DIGEST_LENGTH; ++i){
			if(compare[i] != this->crc[i])
				return false;
		}
		return true;
	}

private:
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry){
		buffer << entry.controller;
		buffer += entry.stride;
		buffer += entry.offset;
		buffer += entry.cLength;
		buffer += entry.uLength;
		buffer += entry.eLength;
		for(U32 i = 0; i < MD5_DIGEST_LENGTH; ++i) buffer += entry.crc[i];
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
		stream.write(reinterpret_cast<const char*>(&entry.crc[0]),    sizeof(uint8_t)*MD5_DIGEST_LENGTH);
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
		for(U32 i = 0; i < MD5_DIGEST_LENGTH; ++i) buffer >> entry.crc[i];
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
		stream.read(reinterpret_cast<char*>(&entry.crc[0]),     sizeof(uint8_t)*MD5_DIGEST_LENGTH);
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
	uint8_t crc[MD5_DIGEST_LENGTH];  // MD5 checksum
	S32 global_key;             // global key
};

}
}



#endif /* CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_ */
