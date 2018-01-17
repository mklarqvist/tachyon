#include <cassert>

#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"
#include "../third_party/zlib/zconf.h"
#include "../third_party/zlib/zlib.h"
#include "../io/BasicBuffer.h"
#include "../algorithm/OpenHashTable.h"
#include "StreamContainer.h"

namespace Tachyon{
namespace Core{

StreamContainer::StreamContainer() :
	n_entries(0),
	n_additions(0)
{}

StreamContainer::StreamContainer(const U32 start_size) :
	n_entries(0),
	n_additions(0),
	buffer_data(start_size),
	buffer_strides(start_size),
	buffer_data_uncompressed(start_size),
	buffer_strides_uncompressed(start_size)
{}

StreamContainer::~StreamContainer(){
	this->buffer_data.deleteAll();
	this->buffer_strides.deleteAll();
	this->buffer_data_uncompressed.deleteAll();
	this->buffer_strides_uncompressed.deleteAll();
}

void StreamContainer::reset(void){
	this->n_entries   = 0;
	this->n_additions = 0;
	this->buffer_data.reset();
	this->buffer_data_uncompressed.reset();
	this->buffer_strides.reset();
	this->buffer_strides_uncompressed.reset();
	this->header.reset();
	this->header_stride.reset();
	this->buffer_data_uncompressed.reset();
	this->buffer_strides_uncompressed.reset();
}

void StreamContainer::resize(const U32 size){
	this->buffer_data.resize(size);
	this->buffer_data_uncompressed.resize(size);
	this->buffer_strides.resize(size);
	this->buffer_strides_uncompressed.resize(size);
}

const bool StreamContainer::generateCRC(void){
	if(this->buffer_data_uncompressed.size() == 0){
		this->header.crc = 0;
	} else {
		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data_uncompressed.data, this->buffer_data_uncompressed.pointer);
		this->header.crc = crc;
	}

	if(this->header.controller.mixedStride == true){
		if(this->buffer_data_uncompressed.size() == 0){
			this->header_stride.crc = 0;
		} else {
			// Checksum for strides
			U32 crc = crc32(0, NULL, 0);
			if(this->buffer_strides_uncompressed.size() > 0){
				crc = crc32(0, NULL, 0);
				crc = crc32(crc, (Bytef*)this->buffer_strides_uncompressed.data, this->buffer_strides_uncompressed.pointer);
				this->header_stride.crc = crc;
			}
		}
	}
	return true;
}

bool StreamContainer::checkCRC(int target){
	if(target == 0){
		if(this->buffer_data_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data_uncompressed.data, this->buffer_data_uncompressed.pointer);
		return(crc == this->header.crc);
	} else if(target == 1){
		if(this->buffer_strides_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_strides_uncompressed.data, this->buffer_strides_uncompressed.pointer);
		return(crc == this->header_stride.crc);
	} else if(target == 3){
		if(this->buffer_data.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data.data, this->buffer_data.pointer);
		return(crc == this->header.crc);
	} else if(target == 4){
		if(this->buffer_strides.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_strides.data, this->buffer_strides.pointer);
		return(crc == this->header_stride.crc);
	}
	return true;
}

bool StreamContainer::checkUniformity(void){
	if(this->n_entries == 0)
		return false;

	// We know the data cannot be uniform if
	// the stride size is uneven
	const S16& stride_size = this->header.stride;
	if(stride_size == -1)
		return false;

	U32 stride_update = stride_size;

	BYTE word_width = sizeof(char);
	switch(this->header.controller.type){
	case YON_TYPE_32B:   stride_update *= sizeof(S32);   word_width = sizeof(S32);   break;
	case YON_TYPE_FLOAT: stride_update *= sizeof(float); word_width = sizeof(float); break;
	case YON_TYPE_8B:    stride_update *= sizeof(char);  word_width = sizeof(char);  break;
	case YON_TYPE_CHAR:  stride_update *= sizeof(char);  word_width = sizeof(char);  break;
	default: return false; break;
	}

	const U64 first_hash = XXH64(&this->buffer_data_uncompressed.data[0], stride_update, 2147483647);

	for(U32 i = 1; i < this->n_entries; ++i){
		if(XXH64(&this->buffer_data_uncompressed.data[i*stride_update], stride_update, 2147483647) != first_hash){
			//std::cerr << "not uniform" << std::endl;
			return(false);
		}
	}

	this->n_entries = 1;
	this->n_additions = 1;
	// Data pointers are updated in case there is no reformatting
	// see StreamContainer::reformat()
	this->buffer_data_uncompressed.pointer = stride_size * word_width;
	this->header.uLength = stride_size * word_width;
	this->header.cLength = stride_size * word_width;
	this->header.controller.uniform = true;
	this->header.controller.mixedStride = false;
	this->header.controller.encoder = Core::YON_ENCODE_NONE;
	return(true);
}

/*
 This function is called during import to shrink each
 word-type to fit min(x) and max(x) in the worst case.
 At this stage all integer values in the stream is of
 type S32. No other values can be shrunk
 */
void StreamContainer::reformat(){
	if(this->buffer_data_uncompressed.pointer)
		return;

	// Recode integer types only
	if(!(this->header.controller.type == YON_TYPE_32B && this->header.controller.signedness == 1))
		return;

	if(this->header.controller.uniform == true)
		return;

	// At this point all integers are S32
	const S32* const dat  = reinterpret_cast<const S32* const>(this->buffer_data_uncompressed.data);
	const U32* const udat = reinterpret_cast<const U32* const>(this->buffer_data_uncompressed.data);
	S32 min = dat[0];
	S32 max = dat[0];
	bool hasMissing = false;

	for(U32 j = 0; j < this->n_additions; ++j){
		if(udat[j] == 0x80000000 || udat[j] == 0x80000001){
			hasMissing = true;
			continue;
		}
		if(dat[j] < min) min = dat[j];
		if(dat[j] > max) max = dat[j];
	}

	BYTE byte_width = 0;
	// If we have missing values then we have to
	// use signedness
	if(min < 0 || hasMissing == true){
		byte_width = ceil((ceil(log2(abs(min) + 1 + 2)) + 1) / 8);  // One bit is used for sign, 2 values for missing and end-of-vector
		const BYTE byte_width_max = ceil((ceil(log2(abs(max) + 1 + 2)) + 1) / 8);
		if(byte_width_max > byte_width){
			byte_width = byte_width_max;
		}
	}
	else byte_width = ceil(ceil(log2(max + 1)) / 8);

	if(byte_width == 0) byte_width = 1;
	else if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if(byte_width > 4) byte_width = 8;

	//std::cerr << "Reformat\t" << min << '\t' << max << '\t' << (int)byte_width << '\t' << (int)hasMissing << std::endl;

	// Phase 2
	// Here we re-encode values using the smallest possible
	// word-size
	this->buffer_data.reset();
	this->buffer_data.resize(this->buffer_data_uncompressed);

	// Is non-negative
	// Also cannot have missing values
	if(min >= 0 && hasMissing == false){
		this->header.controller.signedness = 0;

		if(byte_width == 1){
			this->header.controller.type = YON_TYPE_8B;

			for(U32 j = 0; j < this->n_additions; ++j)
				this->buffer_data += (BYTE)dat[j];

		} else if(byte_width == 2){
			this->header.controller.type = YON_TYPE_16B;

			for(U32 j = 0; j < this->n_additions; ++j)
				this->buffer_data += (U16)dat[j];

		} else if(byte_width == 4){
			this->header.controller.type = YON_TYPE_32B;

			for(U32 j = 0; j < this->n_additions; ++j)
				this->buffer_data += (U32)dat[j];

		} else if(byte_width == 8){
			this->header.controller.type = YON_TYPE_64B;

			for(U32 j = 0; j < this->n_additions; ++j)
				this->buffer_data += (U64)dat[j];

		} else {
			std::cerr << "illegal" << std::endl;
			exit(1);
		}
	}
	// Is negative
	else {
		this->header.controller.signedness = 1;

		if(byte_width == 1){
			this->header.controller.type = YON_TYPE_8B;

			for(U32 j = 0; j < this->n_additions; ++j){
				if(udat[j] == 0x80000000) this->buffer_data += (BYTE)0x80;
				else if(udat[j] == 0x80000001) this->buffer_data += (BYTE)0x81;
				else this->buffer_data += (SBYTE)dat[j];
			}

		} else if(byte_width == 2){
			this->header.controller.type = YON_TYPE_16B;

			for(U32 j = 0; j < this->n_additions; ++j){
				//this->buffer_data += (S16)dat[j];
				if(udat[j] == 0x80000000) this->buffer_data += (U16)0x8000;
				else if(udat[j] == 0x80000001) this->buffer_data += (U16)0x8001;
				else this->buffer_data += (S16)dat[j];
			}

		} else if(byte_width == 4){
			this->header.controller.type = YON_TYPE_32B;

			for(U32 j = 0; j < this->n_additions; ++j){
				//this->buffer_data += (S32)dat[j];
				if(udat[j] == 0x80000000) this->buffer_data += (U32)0x80000000;
				else if(udat[j] == 0x80000001) this->buffer_data += (U32)0x80000001;
				else this->buffer_data += (S32)dat[j];
			}

		} else {
			std::cerr << "illegal" << std::endl;
			exit(1);
		}
	}

	memcpy(this->buffer_data_uncompressed.data, this->buffer_data.data, this->buffer_data.pointer);
	this->buffer_data_uncompressed.pointer = this->buffer_data.pointer;
	this->header.uLength = this->buffer_data_uncompressed.pointer;
	this->buffer_data.reset();
}

void StreamContainer::reformatStride(){
	if(this->buffer_strides_uncompressed.pointer == 0)
		return;

	if(this->header.controller.mixedStride == false)
		return;

	// Recode integer types
	if(!(this->header_stride.controller.type == YON_TYPE_32B &&
	   this->header_stride.controller.signedness == 0)){
		std::cerr << "illegal at this point: " << this->header_stride.controller.type << ":" << this->header_stride.controller.signedness << std::endl;
		exit(1);
	}

	// At this point all integers are S32
	const U32* const dat = reinterpret_cast<const U32* const>(this->buffer_strides_uncompressed.data);
	U32 max = dat[0];

	for(U32 j = 1; j < this->n_entries; ++j)
		if(dat[j] > max) max = dat[j];

	BYTE byte_width = ceil(ceil(log2(max + 1))/8);

	if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if(byte_width > 4) byte_width = 8;
	if(byte_width == 0) byte_width = 1;

	// This cannot ever be uniform
	this->buffer_strides.reset();
	this->buffer_strides.resize(this->buffer_strides_uncompressed.pointer);

	if(byte_width == 1){
		this->header_stride.controller.type = YON_TYPE_8B;

		for(U32 j = 0; j < this->n_entries; ++j)
			this->buffer_strides += (BYTE)dat[j];

	} else if(byte_width == 2){
		this->header_stride.controller.type = YON_TYPE_16B;

		for(U32 j = 0; j < this->n_entries; ++j)
			this->buffer_strides += (U16)dat[j];

	} else if(byte_width == 4){
		this->header_stride.controller.type = YON_TYPE_32B;

		for(U32 j = 0; j < this->n_entries; ++j)
			this->buffer_strides += (U32)dat[j];

	} else if(byte_width == 8){
		this->header_stride.controller.type = YON_TYPE_64B;

		for(U32 j = 0; j < this->n_entries; ++j)
			this->buffer_strides += (U64)dat[j];

	} else {
		std::cerr << "illegal" << std::endl;
		exit(1);
	}
	//std::cerr << "recode shrink strides: " << this->buffer_strides.pointer << '\t' << buffer.pointer << std::endl;
	memcpy(this->buffer_strides_uncompressed.data, this->buffer_strides.data, this->buffer_strides.pointer);
	this->buffer_strides_uncompressed.pointer = this->buffer_strides.pointer;
	this->header_stride.uLength = this->buffer_strides_uncompressed.pointer;
	this->buffer_strides.reset();
}

}
}
