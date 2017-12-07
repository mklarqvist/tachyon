#ifndef CORE_STREAMCONTAINER_H_
#define CORE_STREAMCONTAINER_H_

#include "../io/BasicBuffer.h"
#include "base/StreamContainerHeaderController.h"
#include "base/StreamContainerHeader.h"

namespace Tachyon{
namespace Core{

enum CORE_TYPE{TYPE_8B,
	           TYPE_16B,
			   TYPE_32B,
			   TYPE_64B,
			   TYPE_FLOAT,
			   TYPE_DOUBLE,
			   TYPE_BOOLEAN,
			   TYPE_CHAR,
			   TYPE_STRUCT};

enum CORE_COMPRESSION{ENCODE_NONE,
	                  ENCODE_ZSTD,
					  ENCODE_FSE,
					  ENCODE_RANS};


// Stream container for importing
class StreamContainer{
	typedef StreamContainer self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef StreamContainerHeader header_type;
	typedef StreamContainerHeaderStride header_stride_type;

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

	~StreamContainer(){
		this->buffer_data.deleteAll();
		this->buffer_strides.deleteAll();
		this->buffer_data_uncompressed.deleteAll();
		this->buffer_strides_uncompressed.deleteAll();
	}

	inline void setType(const CORE_TYPE value){ this->header.controller.type = value; }
	inline void setStrideSize(const S32 value){ this->header.stride = value; }
	inline const bool checkStrideSize(const S32& value) const{ return this->header.stride == value; }
	inline void setMixedStrides(void){ this->header.stride = -1; this->header.controller.mixedStride = true; }

	inline void operator++(void){ ++this->n_entries; }

	inline void addStride(const U32& value){ this->buffer_strides += (U32)value; }

	inline void operator+=(const SBYTE& value){
		//assert(this->header.controller.type == 0);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const BYTE& value){
		//assert(this->header.controller.type == 1);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const S16& value){
		//assert(this->header.controller.type == 2);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const U16& value){
		//assert(this->header.controller.type == 3);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const S32& value){
		//assert(this->header.controller.type == 4);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const U32& value){
		//assert(this->header.controller.type == 5);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const U64& value){
		//assert(this->header.controller.type == 6);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const float& value){
		//assert(this->header.controller.type == 7);
		this->buffer_data += value;
		++this->n_additions;
	}

	inline void operator+=(const double& value){
		//assert(this->header.controller.type == 8);
		this->buffer_data += value;
		++this->n_additions;
	}

	void reset(void){
		this->n_entries = 0;
		this->n_additions = 0;
		this->buffer_data.reset();
		this->buffer_strides.reset();
		this->header.reset();
		this->header_stride.reset();
		this->buffer_data_uncompressed.reset();
		this->buffer_strides_uncompressed.reset();
	}

	inline void resize(const U32 size){
		this->buffer_data.resize(size);
		this->buffer_strides.resize(size);
	}

	bool generateCRC(bool both = false);
	bool checkCRC(int target = 0);
	bool checkUniformity(void);
	void reformat(buffer_type& buffer);
	void reformatStride(buffer_type& buffer);
	bool read(std::ifstream& stream, buffer_type& temp_buffer);

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.header;
		if(entry.header.controller.mixedStride)
			stream << entry.header_stride;

		stream << entry.buffer_data;
		if(entry.header.controller.mixedStride)
			stream << entry.buffer_strides;

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.header;
		if(entry.header.controller.mixedStride)
			stream >> entry.header_stride;

		entry.buffer_data.resize(entry.header.uLength);
		stream.read(entry.buffer_data.data, entry.header.cLength);
		entry.buffer_data.pointer = entry.header.cLength;
		if(entry.header.controller.mixedStride){
			entry.buffer_strides.resize(entry.header_stride.uLength);
			stream.read(entry.buffer_strides.data, entry.header_stride.cLength);
			entry.buffer_strides.pointer = entry.header_stride.cLength;
		}
		return(stream);
	}

	inline const U32 getDiskSize(void) const{
		U32 total_size = 0;
		total_size += header.getDiskSize();
		if(this->header.controller.mixedStride)
			total_size += header_stride.getDiskSize();

		total_size += this->buffer_data.pointer;
		if(this->header.controller.mixedStride)
			total_size += this->buffer_strides.pointer;

		return(total_size);
	}

public:
	header_type header;
	header_stride_type header_stride;
	// Not written - used internally only during import
	U32 n_entries;   // number of container entries
	U32 n_additions; // number of times an addition operation was executed
	// Buffers - only bit that are written to disk
	// from here
	buffer_type buffer_data;
	buffer_type buffer_strides;
	// These buffers are for internal use only
	// They are used during decompression and are
	// not written to disk
	buffer_type buffer_data_uncompressed;
	buffer_type buffer_strides_uncompressed;
};

}
}

#endif /* CORE_STREAMCONTAINER_H_ */
