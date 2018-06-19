#ifndef CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_
#define CONTAINERS_COMPONENTS_DATA_CONTAINER_HEADER_OBJECT_H_

namespace tachyon{
namespace containers{

struct DataContainerHeaderObject{
	typedef DataContainerHeaderObject     self_type;
	typedef DataContainerHeaderController controller_type;

	DataContainerHeaderObject() :
		stride(1),
		offset(0),
		cLength(0),
		uLength(0),
		eLength(0),
		crc(0),
		global_key(-1)
	{}

	DataContainerHeaderObject(const DataContainerHeaderObject& other) :
		controller(other.controller),
		stride(other.stride),
		offset(other.offset),
		cLength(other.cLength),
		uLength(other.uLength),
		eLength(other.eLength),
		crc(other.crc),
		global_key(other.global_key)
	{
	}

	/* noexcept needed to enable optimizations in containers */
	DataContainerHeaderObject(DataContainerHeaderObject&& other) noexcept :
		controller(other.controller),
		stride(other.stride),
		offset(other.offset),
		cLength(other.cLength),
		uLength(other.uLength),
		eLength(other.eLength),
		crc(other.crc),
		global_key(other.global_key)
	{

	}

	 // copy assignment
	DataContainerHeaderObject& operator=(const DataContainerHeaderObject& other){
		this->controller = other.controller;
		this->stride     = other.stride;
		this->offset     = other.offset;
		this->cLength    = other.cLength;
		this->uLength    = other.uLength;
		this->eLength    = other.eLength;
		this->crc        = other.crc;
		this->global_key = other.global_key;
		return *this;
	}


	/** Move assignment operator */
	DataContainerHeaderObject& operator=(DataContainerHeaderObject&& other) noexcept{
		this->controller = other.controller;
		this->stride     = other.stride;
		this->offset     = other.offset;
		this->cLength    = other.cLength;
		this->uLength    = other.uLength;
		this->eLength    = other.eLength;
		this->crc        = other.crc;
		this->global_key = other.global_key;
		return *this;
	}

	~DataContainerHeaderObject(){ }

	inline void reset(void){
		this->controller.clear();
		this->stride     = 1;
		this->offset     = 0;
		this->cLength    = 0;
		this->uLength    = 0;
		this->crc        = 0;
		this->global_key = -1;
	}

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

	const bool operator==(const self_type& other) const{
		if(this->stride     != other.stride)     return false;
		if(this->offset     != other.offset)     return false;
		if(this->cLength    != other.cLength)    return false;
		if(this->uLength    != other.uLength)    return false;
		if(this->eLength    != other.eLength)    return false;
		if(this->crc        != other.crc)        return false;
		if(this->global_key != other.global_key) return false;
		if(this->controller != other.controller) return false;
		return true;
	}

	inline const bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	const SBYTE getPrimitiveWidth(void) const{
		// We do not care about signedness here
		switch(this->controller.type){
		case(YON_TYPE_UNKNOWN):
		case(YON_TYPE_STRUCT): return(-1);
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):   return(sizeof(char));
		case(YON_TYPE_8B):     return(sizeof(BYTE));
		case(YON_TYPE_16B):    return(sizeof(U16));
		case(YON_TYPE_32B):    return(sizeof(U32));
		case(YON_TYPE_64B):    return(sizeof(U64));
		case(YON_TYPE_FLOAT):  return(sizeof(float));
		case(YON_TYPE_DOUBLE): return(sizeof(double));
		}
		return 0;
	}

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
