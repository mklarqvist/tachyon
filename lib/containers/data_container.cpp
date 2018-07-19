#include <cassert>

#include "data_container.h"
#include "algorithm/OpenHashTable.h"
#include "io/basic_buffer.h"
#include "support/helpers.h"
#include "support/type_definitions.h"
#include "third_party/zlib/zconf.h"
#include "third_party/zlib/zlib.h"

namespace tachyon{
namespace containers{

DataContainer::DataContainer()
{}

DataContainer::DataContainer(const U32 start_size) :
	buffer_data(start_size),
	buffer_strides(start_size),
	buffer_data_uncompressed(start_size),
	buffer_strides_uncompressed(start_size)
{}

DataContainer::~DataContainer(){ }

void DataContainer::reset(void){
	this->buffer_data.reset();
	this->buffer_data_uncompressed.reset();
	this->buffer_strides.reset();
	this->buffer_strides_uncompressed.reset();
	this->header.reset();
}

void DataContainer::resize(const U32 size){
	this->buffer_data.resize(size);
	this->buffer_data_uncompressed.resize(size);
	this->buffer_strides.resize(size);
	this->buffer_strides_uncompressed.resize(size);
}

void DataContainer::generateCRC(void){
	if(this->buffer_data_uncompressed.size() == 0){
		this->header.data_header.crc = 0;
	} else {
		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data_uncompressed.buffer, this->buffer_data_uncompressed.size());
		this->header.data_header.crc = crc;
	}

	if(this->header.data_header.hasMixedStride()){
		if(this->buffer_data_uncompressed.size() == 0){
			this->header.stride_header.crc = 0;
		} else {
			// Checksum for strides
			U32 crc = crc32(0, NULL, 0);
			if(this->buffer_strides_uncompressed.size() > 0){
				crc = crc32(0, NULL, 0);
				crc = crc32(crc, (Bytef*)this->buffer_strides_uncompressed.buffer, this->buffer_strides_uncompressed.size());
				this->header.stride_header.crc = crc;
			}
		}
	}
}

bool DataContainer::checkCRC(int target){
	if(target == 0){
		if(this->buffer_data_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data_uncompressed.buffer, this->buffer_data_uncompressed.size());
		return(crc == this->header.data_header.crc);
	} else if(target == 1){
		if(this->buffer_strides_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_strides_uncompressed.buffer, this->buffer_strides_uncompressed.size());
		return(crc == this->header.stride_header.crc);
	} else if(target == 3){
		if(this->buffer_data.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data.buffer, this->buffer_data.size());
		return(crc == this->header.data_header.crc);
	} else if(target == 4){
		if(this->buffer_strides.size() == 0)
			return true;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_strides.buffer, this->buffer_strides.size());
		return(crc == this->header.stride_header.crc);
	}
	return true;
}

bool DataContainer::checkUniformity(void){
	if(this->header.n_entries == 0)
		return false;

	// We know the stride cannot be uniform if
	// the stride size is uneven
	const S16& stride_size = this->header.data_header.stride;
	if(stride_size == -1)
		return false;

	U32 stride_update = stride_size;

	BYTE word_width = sizeof(char);
	switch(this->header.data_header.controller.type){
	case YON_TYPE_DOUBLE: stride_update *= sizeof(double); word_width = sizeof(double); break;
	case YON_TYPE_FLOAT:  stride_update *= sizeof(float);  word_width = sizeof(float);  break;
	case YON_TYPE_8B:     stride_update *= sizeof(BYTE);   word_width = sizeof(BYTE);   break;
	case YON_TYPE_16B:    stride_update *= sizeof(U16);    word_width = sizeof(U16);    break;
	case YON_TYPE_32B:    stride_update *= sizeof(S32);    word_width = sizeof(S32);    break;
	case YON_TYPE_64B:    stride_update *= sizeof(U64);    word_width = sizeof(U64);    break;
	case YON_TYPE_CHAR:   stride_update *= sizeof(char);   word_width = sizeof(char);   break;
	default: return false; break;
	}

	const U64 first_hash = XXH64(this->buffer_data_uncompressed.data(), stride_update, 2147483647);

	U64 cumulative_position = stride_update;
	for(U32 i = 1; i < this->header.n_entries; ++i){
		if(XXH64(&this->buffer_data_uncompressed.buffer[cumulative_position], stride_update, 2147483647) != first_hash){
			return(false);
		}
		cumulative_position += stride_update;
	}
	assert(cumulative_position == this->buffer_data_uncompressed.size());

	this->header.n_entries   = 1;
	this->header.n_strides   = 0;
	this->buffer_data_uncompressed.n_chars          = stride_size * word_width;
	this->header.data_header.uLength                = stride_size * word_width;
	this->header.data_header.cLength                = stride_size * word_width;
	this->header.data_header.controller.uniform     = true;
	this->header.data_header.controller.mixedStride = false;
	this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
	return(true);
}

void DataContainer::reformatInteger(){
	if(this->buffer_data_uncompressed.size() == 0)
		return;

	// Recode integer types only
	if(!(this->header.data_header.controller.type == YON_TYPE_32B && this->header.data_header.controller.signedness == 1)){
		return;
	}

	// Do not recode if the data is uniform
	if(this->header.data_header.isUniform())
		return;

	// At this point all integers are S32
	const S32* const dat  = reinterpret_cast<const S32* const>(this->buffer_data_uncompressed.data());
	const U32* const udat = reinterpret_cast<const U32* const>(this->buffer_data_uncompressed.data());
	S32 min = dat[0];
	S32 max = dat[0];
	bool hasMissing = false;

	for(U32 j = 0; j < this->header.n_additions; ++j){
		if(udat[j] == 0x80000000 || udat[j] == 0x80000001){
			hasMissing = true;
			continue;
		}
		if(dat[j] < min) min = dat[j];
		if(dat[j] > max) max = dat[j];
	}

	BYTE byte_width = 0;
	// If we have missing values then we have to use signedness
	if(min < 0 || hasMissing == true){
		byte_width = ceil((ceil(log2(abs(min) + 1 + 2)) + 1) / 8);  // One bit is used for sign, 2 values for missing and end-of-vector
		const BYTE byte_width_max = ceil((ceil(log2(abs(max) + 1 + 2)) + 1) / 8);
		if(byte_width_max > byte_width){
			byte_width = byte_width_max;
		}
	}
	else byte_width = ceil(ceil(log2(max + 1)) / 8);

	if(byte_width == 0)     byte_width = 1;
	else if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if(byte_width > 4) byte_width = 8;

	// Phase 2
	// Here we re-encode values using the smallest possible word-size
	this->buffer_data.reset();
	this->buffer_data.resize(this->buffer_data_uncompressed.size() + 65536);

	// Is non-negative
	// Also cannot have missing values
	if(min >= 0 && hasMissing == false){
		this->header.data_header.controller.signedness = 0;

		if(byte_width == 1){
			this->header.data_header.controller.type = YON_TYPE_8B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				this->buffer_data += (BYTE)dat[j];
			}

		} else if(byte_width == 2){
			this->header.data_header.controller.type = YON_TYPE_16B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				this->buffer_data += (U16)dat[j];
			}

		} else if(byte_width == 4){
			this->header.data_header.controller.type = YON_TYPE_32B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				this->buffer_data += (U32)dat[j];
			}

		} else if(byte_width == 8){
			this->header.data_header.controller.type = YON_TYPE_64B;

			for(U32 j = 0; j < this->header.n_additions; ++j)
				this->buffer_data += (U64)dat[j];

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}
	// Is negative or has missing
	else {
		this->header.data_header.controller.signedness = true;

		if(byte_width == 1){
			this->header.data_header.controller.type = YON_TYPE_8B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(udat[j] == 0x80000000)      this->buffer_data += (BYTE)0x80;
				else if(udat[j] == 0x80000001) this->buffer_data += (BYTE)0x81;
				else this->buffer_data += (SBYTE)dat[j];
			}

		} else if(byte_width == 2){
			this->header.data_header.controller.type = YON_TYPE_16B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(udat[j] == 0x80000000)      this->buffer_data += (U16)0x8000;
				else if(udat[j] == 0x80000001) this->buffer_data += (U16)0x8001;
				else this->buffer_data += (S16)dat[j];
			}

		} else if(byte_width == 4){
			this->header.data_header.controller.type = YON_TYPE_32B;

			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(udat[j] == 0x80000000)      this->buffer_data += (U32)0x80000000;
				else if(udat[j] == 0x80000001) this->buffer_data += (U32)0x80000001;
				else this->buffer_data += (S32)dat[j];
			}

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}

	memcpy(this->buffer_data_uncompressed.buffer, this->buffer_data.buffer, this->buffer_data.size());
	this->buffer_data_uncompressed.n_chars = this->buffer_data.size();
	this->header.data_header.uLength       = this->buffer_data_uncompressed.size();
	this->buffer_data.reset();
}

void DataContainer::reformatStride(){
	if(this->buffer_strides_uncompressed.size() == 0)
		return;

	if(this->header.data_header.hasMixedStride() == false)
		return;

	// Recode integer types
	if(!(this->header.stride_header.controller.type == YON_TYPE_32B &&
	   this->header.stride_header.controller.signedness == 0)){
		std::cerr << utility::timestamp("ERROR") << "Illegal to have non-U32 values at this point: " << this->header.stride_header.controller.type << ":" << this->header.stride_header.controller.signedness << std::endl;
		exit(1);
	}

	// At this point all integers are U32
	const U32* const dat = reinterpret_cast<const U32* const>(this->buffer_strides_uncompressed.data());
	assert(this->buffer_strides_uncompressed.size() % sizeof(U32) == 0);
	assert(this->buffer_strides_uncompressed.size() / sizeof(U32) == this->header.n_strides);

	U32 max = 1;
	for(U32 j = 0; j < this->header.n_strides; ++j){
		if(dat[j] > max)
			max = dat[j];
	}

	BYTE byte_width = ceil(log2(max + 1) / 8);
	if(byte_width == 0) byte_width = 1;
	if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	if(byte_width > 4)  byte_width = 8;

	// This cannot ever be uniform
	this->buffer_strides.reset();
	this->buffer_strides.resize(this->buffer_strides_uncompressed.size() + 65536);

	if(byte_width == 1){
		this->header.stride_header.controller.type = YON_TYPE_8B;

		for(U32 j = 0; j < this->header.n_strides; ++j){
			//assert((BYTE)dat[j] == dat[j]);
			this->buffer_strides += (BYTE)dat[j];
		}

	} else if(byte_width == 2){
		this->header.stride_header.controller.type = YON_TYPE_16B;

		for(U32 j = 0; j < this->header.n_strides; ++j){
			//assert((U16)dat[j] == dat[j]);
			this->buffer_strides += (U16)dat[j];
		}

	} else if(byte_width == 4){
		this->header.stride_header.controller.type = YON_TYPE_32B;

		for(U32 j = 0; j < this->header.n_strides; ++j){
			//assert((U32)dat[j] == dat[j]);
			this->buffer_strides += (U32)dat[j];
		}

	} else if(byte_width == 8){
		this->header.stride_header.controller.type = YON_TYPE_64B;

		for(U32 j = 0; j < this->header.n_strides; ++j){
			//assert((U64)dat[j] == dat[j]);
			this->buffer_strides += (U64)dat[j];
		}

	} else {
		std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
		exit(1);
	}
	memcpy(this->buffer_strides_uncompressed.data(), this->buffer_strides.data(), this->buffer_strides.size());
	this->buffer_strides_uncompressed.n_chars = this->buffer_strides.size();
	this->header.stride_header.uLength        = this->buffer_strides_uncompressed.size();
	this->buffer_strides.reset();
}

U32 DataContainer::getObjectSize(void) const{
	// In case data is encrypted
	if(this->header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
		return(this->buffer_data.size());

	U32 total_size = this->buffer_data.size();
	if(this->header.data_header.hasMixedStride())
		total_size += this->buffer_strides.size();

	return(total_size);
}

U64 DataContainer::getObjectSizeUncompressed(void) const{
	U64 total_size = this->buffer_data_uncompressed.size();
	if(this->header.data_header.hasMixedStride())
		total_size += this->buffer_strides_uncompressed.size();

	return(total_size);
}

void DataContainer::deltaEncode(){
	if(this->size() == 0)
		return;

	// Recode integer types only
	if(!(this->header.data_header.controller.type == YON_TYPE_32B && this->header.data_header.controller.signedness == 1)){
		return;
	}

	if(this->header.data_header.controller.uniform == true)
		return;

	// At this point all integers are S32
	const S32* const dat  = reinterpret_cast<const S32* const>(this->buffer_data_uncompressed.buffer);

	// check for uniformity except first
	if(this->header.n_additions > 1){
		bool is_uniform_delta = true;
		const S32 test_diff = dat[1] - dat[0];
		for(U32 i = 2; i < this->header.n_additions; ++i){
			if(dat[i] - dat[i - 1] != test_diff){
				is_uniform_delta = false;
				break;
			}
		}

		if(is_uniform_delta){
			this->header.n_entries   = 1;
			this->header.n_additions = 1;
			// Data pointers are updated in case there is no reformatting
			// see StreamContainer::reformat()
			this->buffer_data_uncompressed.n_chars             = sizeof(S32);
			this->header.data_header.uLength                   = sizeof(S32);
			this->header.data_header.cLength                   = sizeof(S32);
			this->header.data_header.controller.uniform        = true;
			this->header.data_header.controller.mixedStride    = false;
			this->header.data_header.controller.encoder        = YON_ENCODE_NONE;
			return;
		}
	}

	this->buffer_data += dat[0];
	for(U32 j = 1; j < this->header.n_additions; ++j){
		this->buffer_data += dat[j] - dat[j-1];
	}
	memcpy(this->buffer_data_uncompressed.data(),
			this->buffer_data.data(),
			this->buffer_data.size());

}

void DataContainer::updateContainer(bool reformat){
	// If the data container has entries in it but has
	// no actual data then it is a BOOLEAN
	if(this->header.n_entries && this->buffer_data_uncompressed.size() == 0){
		this->header.reset();
		this->header.data_header.controller.type        = YON_TYPE_BOOLEAN;
		this->header.data_header.controller.uniform     = true;
		this->header.data_header.controller.mixedStride = false;
		this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
		this->header.data_header.controller.signedness  = 0;
		this->header.data_header.stride  = 0;
		this->header.data_header.uLength = 0;
		this->header.data_header.cLength = 0;
		this->header.n_strides           = 0;
		return;
	}

	if(this->buffer_data_uncompressed.size() == 0)
		return;

	// Check if stream is uniform in content
	if(this->header.data_header.controller.type != YON_TYPE_STRUCT){
		this->checkUniformity();
		// Reformat stream to use as small word size as possible
		if(reformat) this->reformatInteger();
	}

	// Set uncompressed length
	this->header.data_header.uLength = this->buffer_data_uncompressed.size();

	// If we have mixed striding
	if(this->header.data_header.hasMixedStride()){
		// Reformat stream to use as small word size as possible
		if(reformat) this->reformatStride();
		this->header.stride_header.uLength = this->buffer_strides_uncompressed.size();
	}
}

}
}
