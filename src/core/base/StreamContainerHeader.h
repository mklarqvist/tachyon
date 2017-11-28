#ifndef CORE_BASE_STREAMCONTAINERHEADER_H_
#define CORE_BASE_STREAMCONTAINERHEADER_H_

#include <iostream>
#include <fstream>
#include <cstring>

#include "../../support/TypeDefinitions.h"

namespace Tachyon{
namespace Core{

/*
 Stream header data structure
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+----+
 |   2   |   2   |       4       |        4      |       4       |       4       | X ~
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 ^       ^       ^               ^               ^               ^               ^
 |       |       |               |               |               |               |
 CNT     STRIDE  |               COMP LENGTH     UNCOMP LENGTH   CRC             POSSIBLE COMPRESSION PARAMETERS
                 OFFSET
*/
struct StreamContainerHeader{
	typedef StreamContainerHeader self_type;
	typedef StreamContainerHeaderController controller_type;

	StreamContainerHeader() :
		stride(-1),
		offset(0),
		cLength(0),
		uLength(0),
		crc(0),
		n_extra(0),
		extra(nullptr)
	{}

	StreamContainerHeader(const StreamContainerHeader& other) :
		controller(other.controller),
		stride(other.stride),
		offset(other.offset),
		cLength(other.cLength),
		uLength(other.uLength),
		crc(other.crc),
		n_extra(other.n_extra),
		extra(nullptr)
	{
		if(other.extra != nullptr){
			this->extra = new char[n_extra];
			memcpy(this->extra, other.extra, other.n_extra);
		}
	}

	/* noexcept needed to enable optimizations in containers */
	StreamContainerHeader(StreamContainerHeader&& other) noexcept :
		controller(other.controller),
		stride(other.stride),
		offset(other.offset),
		cLength(other.cLength),
		uLength(other.uLength),
		crc(other.crc),
		n_extra(other.n_extra),
		extra(other.extra)
	{
		other.extra = nullptr;
	}

	 // copy assignment
	StreamContainerHeader& operator=(const StreamContainerHeader& other){
		this->controller = other.controller;
		this->stride = other.stride;
		this->offset = other.offset;
		this->cLength = other.cLength;
		this->uLength = other.uLength;
		this->crc = other.crc;

		if(other.n_extra > 0){
			char* tmp = new char[other.n_extra];
			memcpy(tmp, other.extra, other.n_extra);
			delete[] this->extra;
			this->extra = tmp;
			this->n_extra = other.n_extra;
		} else this->extra = nullptr;
		return *this;
	}


	/** Move assignment operator */
	StreamContainerHeader& operator=(StreamContainerHeader&& other) noexcept{
		this->controller = other.controller;
		this->stride = other.stride;
		this->offset = other.offset;
		this->cLength = other.cLength;
		this->uLength = other.uLength;
		this->crc = other.crc;
		this->n_extra = other.n_extra;
		delete [] this->extra;
		this->extra = other.extra;
		other.extra = nullptr;
		return *this;
	}

	~StreamContainerHeader(){delete [] this->extra;}

	inline void reset(void){
		this->controller.clear();
		this->stride = -1;
		this->offset = 0;
		this->cLength = 0;
		this->uLength = 0;
		this->crc = 0;
		this->n_extra = 0;
		delete [] this->extra;
		this->extra = nullptr;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.stride), sizeof(S16));
		stream.write(reinterpret_cast<const char*>(&entry.offset), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.cLength),sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uLength),sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.crc),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_extra),sizeof(U16));
		if(entry.n_extra > 0){
			stream.write(entry.extra, entry.n_extra);
		}

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.stride),  sizeof(S16));
		stream.read(reinterpret_cast<char*>(&entry.offset),  sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.cLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.crc),     sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_extra), sizeof(U16));
		if(entry.n_extra > 0){
			delete [] entry.extra;
			entry.extra = new char[entry.n_extra];
			stream.read(entry.extra, entry.n_extra);
		}
		return(stream);
	}

	const U32 getDiskSize(void) const{
		U32 total_size = 3*sizeof(U16) + sizeof(U32)*4;
		if(this->n_extra > 0)
			total_size += this->n_extra;

		return total_size;
	}

public:
	controller_type controller;
	S16 stride;
	U32 offset;
	U32 cLength;
	U32 uLength;
	U32 crc;
	U16 n_extra;
	char* extra; // extra length is encoder specific
};

/*
 If stride is mixed then use this
 secondary structure
 +---+---+---+---+---+---+---+---+---+---+---+
 |   2   |       4       |       4       | X ~
 +---+---+---+---+---+---+---+---+---+---+---+
 ^       ^               ^               ^
 |       |               |               |
 CNT     COMP LENGTH     UNCOMP LENGTH   POSSIBLE COMPRESSION PARAMETERS

*/
struct StreamContainerHeaderStride{
	typedef StreamContainerHeaderStride self_type;
	typedef StreamContainerHeaderController controller_type;

	StreamContainerHeaderStride() :
		cLength(0),
		uLength(0),
		crc(0),
		n_extra(0),
		extra(nullptr)
	{}

	StreamContainerHeaderStride(const self_type& other) :
		controller(other.controller),
		cLength(other.cLength),
		uLength(other.uLength),
		crc(other.crc),
		n_extra(other.n_extra),
		extra(nullptr)
	{
		if(other.extra != nullptr){
			this->extra = new char[n_extra];
			memcpy(this->extra, other.extra, other.n_extra);
		}
	}

	/* noexcept needed to enable optimizations in containers */
	StreamContainerHeaderStride(self_type&& other) noexcept :
		controller(other.controller),
		cLength(other.cLength),
		uLength(other.uLength),
		crc(other.crc),
		n_extra(other.n_extra),
		extra(other.extra)
	{
		other.extra = nullptr;
	}

	self_type& operator=(const self_type& other){
		self_type tmp(other); // re-use copy-constructor
		*this = std::move(tmp);                 // re-use move-assignment
		return *this;
	}

	/** Move assignment operator */
	self_type& operator=(self_type&& other) noexcept{
		// prevent self-move
		if(this!=&other){
			this->controller = other.controller;
			this->cLength = other.cLength;
			this->uLength = other.uLength;
			this->crc = other.crc;
			this->n_extra = other.n_extra;

			if(other.extra != nullptr){
				delete [] this->extra;
				this->extra = other.extra;
				other.extra = nullptr;
			}

		}
		return *this;
	}

	~StreamContainerHeaderStride(){
		delete [] this->extra;
	}

	inline void reset(void){
		this->controller.clear();
		this->cLength = 0;
		this->uLength = 0;
		this->crc = 0;
		this->n_extra = 0;
		delete [] this->extra;
		this->extra = nullptr;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.cLength), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uLength), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.crc),     sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_extra), sizeof(U16));
		if(entry.n_extra > 0)
			stream.write(entry.extra, entry.n_extra);

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.cLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.crc),     sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_extra), sizeof(U16));
		if(entry.n_extra > 0){
			delete [] entry.extra;
			entry.extra = new char[entry.n_extra];
			stream.read(entry.extra, entry.n_extra);
		}
		return(stream);
	}

	const U32 getDiskSize(void) const{
		U32 total_size = 2*sizeof(U16) + sizeof(U32)*3;
		if(this->n_extra > 0)
			total_size += this->n_extra;

		return total_size;
	}

public:
	controller_type controller;
	U32 cLength;
	U32 uLength;
	U32 crc;
	U16 n_extra;
	char* extra; // extra length is encoder specific
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADER_H_ */
