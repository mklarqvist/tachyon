#ifndef STREAMCONTAINER_H_
#define STREAMCONTAINER_H_

#include "../io/BasicBuffer.h"

namespace Tomahawk{
namespace Core{

enum CORE_TYPE{TYPE_8B, TYPE_16B, TYPE_32B, TYPE_64B};

struct StreamContainerHeaderController{

	BYTE type: 3,
	     signedness: 1,
		 encoder: 2,
		 other: 2;
};

/*
 Data starts at offset
 +---+---+---+---+---+---+---+---+
 | 1 | 2 | 4 | 4 | 4 | 1 | 4 | 4 |
 +---+---+---+---+---+---+---+---+
 ^   ^   ^   ^   ^   ^   ^   ^
 |   |   |   |   |   |   |   |
 CNT STRIDE  CLENGTH CNT CLENGTH
         OFFSET  ULENGTH     ULENGTH
*/
struct StreamContainerHeader{
	typedef StreamContainerHeader self_type;
	typedef StreamContainerHeaderController controller_type;

	StreamContainerHeader() :
		stride(-1),
		offset(0),
		cLength(0),
		uLength(0)
	{}
	virtual ~StreamContainerHeader(){}

public:
	controller_type controller;
	S16 stride;
	U32 offset;
	U32 cLength;
	U32 uLength;
};

struct StreamContainerHeaderStride : public StreamContainerHeader{
	typedef StreamContainerHeader self_type;

	StreamContainerHeaderStride() :
		controllerStride(0),
		cLengthStride(0),
		uLengthStride(0)
	{}
	~StreamContainerHeaderStride(){}

public:
	BYTE controllerStride;
	U32 cLengthStride;
	U32 uLengthStride;
};

// Stream container for importing
class StreamContainer{
	typedef StreamContainer self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef StreamContainerHeader header_type;

public:
	StreamContainer() :
		n_entries(0),
		n_additions(0)
	{}

	StreamContainer(const U32 start_size) :
		n_entries(0),
		n_additions(0),
		buffer_data(start_size),
		buffer_strides(start_size)
	{}

	~StreamContainer(){ this->buffer_data.deleteAll(); }

	inline void setType(const BYTE value){ this->header.controller.type = value; }
	inline void setStrideSize(const S32 value){ this->header.stride = value; }
	inline const bool checkStrideSize(const S32& value) const{ return this->header.stride == value; }
	inline void setMixedStrides(void){ this->header.stride = -1; }

	inline void operator++(void){ ++this->n_additions; }

	inline void addStride(const U32& value){ this->buffer_strides += value; }

	inline void operator+=(const SBYTE& value){
		if(this->header.controller.type != 0) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const BYTE& value){
		if(this->header.controller.type != 1) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const S16& value){
		if(this->header.controller.type != 2) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U16& value){
		if(this->header.controller.type != 3) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const S32& value){
		if(this->header.controller.type != 4) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U32& value){
		if(this->header.controller.type != 5) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U64& value){
		if(this->header.controller.type != 6) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const float& value){
		if(this->header.controller.type != 7) exit(1);
		this->buffer_data += value;
		++this->n_entries;
	}

	void reset(void){
		this->header.controller.type = 0;
		this->n_entries = 0;
		this->n_additions = 0;
		this->buffer_data.reset();
		this->buffer_strides.reset();
		this->header.stride = -1;
	}

	inline void resize(const U32 size){
		this->buffer_data.resize(size);
		this->buffer_strides.resize(size);
	}

public:
	header_type header;
	U32 n_entries;   // number of entries
	U32 n_additions; // number of times added
	buffer_type buffer_data;
	buffer_type buffer_strides;
};

}
}

#endif /* STREAMCONTAINER_H_ */
