#include <cassert>

#include "data_container.h"
#include "io/basic_buffer.h"
#include "support/helpers.h"
#include "support/type_definitions.h"

#include "algorithm/digest/variant_digest_manager.h"

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

void DataContainer::GenerateMd5(void){
	algorithm::VariantDigestManager::GenerateMd5(this->buffer_data_uncompressed.data(), this->buffer_data_uncompressed.size(), &this->header.data_header.crc[0]);
	algorithm::VariantDigestManager::GenerateMd5(this->buffer_strides_uncompressed.data(), this->buffer_strides_uncompressed.size(), &this->header.stride_header.crc[0]);
}

bool DataContainer::CheckMd5(int target){
	if(target == 0){
		if(this->buffer_data_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->buffer_data_uncompressed.data(), this->buffer_data_uncompressed.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if(target == 1){
		if(this->buffer_strides_uncompressed.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->buffer_strides_uncompressed.data(), this->buffer_strides_uncompressed.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));

	} else if(target == 3){
		if(this->buffer_data.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->buffer_data.data(), this->buffer_data.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if(target == 4){
		if(this->buffer_strides.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->buffer_strides.data(), this->buffer_strides.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));
	}
	return true;
}

bool DataContainer::CheckUniformity(void){
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

void DataContainer::ReformatInteger(){
	if(this->buffer_data_uncompressed.size() == 0)
		return;

	// Recode integer types only.
	if(!(this->header.data_header.controller.type == YON_TYPE_32B &&
	     this->header.data_header.controller.signedness == true))
	{
		return;
	}

	// Do not recode if the data is uniform.
	if(this->header.data_header.IsUniform())
		return;

	// At this point all integers are signed 32-bit integers.
	assert(this->buffer_data_uncompressed.size() % sizeof(S32) == 0);
	assert(this->header.n_additions * sizeof(S32) == this->buffer_data_uncompressed.size());
	const S32* const dat  = reinterpret_cast<const S32* const>(this->buffer_data_uncompressed.data());
	S32 min_value = dat[0];
	S32 max_value = dat[0];
	bool has_special = false;

	// Iterate over available data and search for either missingness
	// or sentinel node values. If a match is found trigger the
	// bool flag to signal the use of the recoding procedure.
	for(U32 j = 0; j < this->header.n_additions; ++j){
		if(dat[j] == bcf_int32_missing || dat[j] == bcf_int32_vector_end){
			has_special = true;
			continue;
		}
		if(dat[j] < min_value) min_value = dat[j];
		if(dat[j] > max_value) max_value = dat[j];
	}

	// If we have missing values then we have to use signed
	// primitives to accommodate this fact.
	BYTE byte_width = 0;
	if(min_value < 0 || has_special == true){
		byte_width = ceil((ceil(log2(abs(min_value) + 1 + 2)) + 1) / 8);  // One bit is used for sign, 2 values for missing and end-of-vector
		const BYTE byte_width_max = ceil((ceil(log2(abs(max_value) + 1 + 2)) + 1) / 8);
		if(byte_width_max > byte_width){
			byte_width = byte_width_max;
		}
	}
	else byte_width = ceil(ceil(log2(max_value + 1)) / 8);

	// Select the smallest primitive type (word width) that
	// can hold the target data range.
	if(byte_width == 0)     byte_width = 1;
	else if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if(byte_width > 4) byte_width = 8;

	// Setup buffers.
	this->buffer_data.reset();
	this->buffer_data.resize(this->buffer_data_uncompressed.size() + 65536);

	// Is non-negative
	// Also cannot have missing values
	if(min_value >= 0 && has_special == false){
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
	// Is negative or Has missing
	else {
		this->header.data_header.controller.signedness = true;

		if(byte_width == 1){
			this->header.data_header.controller.type = YON_TYPE_8B;

			const SBYTE missing = INT8_MIN;
			const SBYTE eov     = INT8_MIN + 1;
			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->buffer_data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->buffer_data += eov;
				else this->buffer_data += (SBYTE)dat[j];
			}

		} else if(byte_width == 2){
			this->header.data_header.controller.type = YON_TYPE_16B;

			const S16 missing = INT16_MIN;
			const S16 eov     = INT16_MIN + 1;
			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->buffer_data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->buffer_data += eov;
				else this->buffer_data += (S16)dat[j];
			}

		} else if(byte_width == 4){
			this->header.data_header.controller.type = YON_TYPE_32B;

			const S32 missing = INT32_MIN;
			const S32 eov     = INT32_MIN + 1;
			for(U32 j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->buffer_data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->buffer_data += eov;
				else this->buffer_data += (S32)dat[j];
			}

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}
	assert(this->buffer_data.size() % byte_width == 0);
	assert(this->header.n_additions * byte_width == this->buffer_data.size());

	memcpy(this->buffer_data_uncompressed.buffer, this->buffer_data.buffer, this->buffer_data.size());
	this->buffer_data_uncompressed.n_chars = this->buffer_data.size();
	this->header.data_header.uLength       = this->buffer_data_uncompressed.size();
	this->buffer_data.reset();
}

void DataContainer::ReformatStride(){
	if(this->buffer_strides_uncompressed.size() == 0)
		return;

	if(this->header.data_header.HasMixedStride() == false)
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

U32 DataContainer::GetObjectSize(void) const{
	// In case data is encrypted
	if(this->header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
		return(this->buffer_data.size());

	U32 total_size = this->buffer_data.size();
	if(this->header.data_header.HasMixedStride())
		total_size += this->buffer_strides.size();

	return(total_size);
}

U64 DataContainer::GetObjectSizeUncompressed(void) const{
	U64 total_size = this->buffer_data_uncompressed.size();
	if(this->header.data_header.HasMixedStride())
		total_size += this->buffer_strides_uncompressed.size();

	return(total_size);
}

void DataContainer::UpdateContainer(bool reformat_data, bool reformat_stride){
	// If the data container Has entries in it but Has
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
		this->CheckUniformity();
		// Reformat stream to use as small word size as possible
		if(reformat_data) this->ReformatInteger();
	}

	// Set uncompressed length
	this->header.data_header.uLength = this->buffer_data_uncompressed.size();

	// If we have mixed striding
	if(this->header.data_header.HasMixedStride()){
		// Reformat stream to use as small word size as possible
		if(reformat_stride) this->ReformatStride();
		this->header.stride_header.uLength = this->buffer_strides_uncompressed.size();
	}
}

void DataContainer::AddStride(const U32 value){
	// If this is the first stride set
	if(this->header.n_strides == 0){
		this->header.stride_header.controller.type = YON_TYPE_32B;
		this->header.stride_header.controller.signedness = false;
		this->SetStrideSize(value);
	}

	// Check if there are different strides
	if(!this->CheckStrideSize(value)){
		this->TriggerMixedStride();
	}

	// Add value
	this->buffer_strides_uncompressed += (U32)value;
	++this->header.n_strides;
}

bool DataContainer::Add(const BYTE& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added BYTE" << std::endl;
		exit(1);
		return false;
	}
	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const U16& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added U16" << std::endl;
		exit(1);
		return false;
	}
	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const U32& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()){
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added U32" << std::endl;
		exit(1);
		return false;
	}
	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const SBYTE& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}

	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added SBYTE" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int8_vector_end){
		this->buffer_data_uncompressed += (S32)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int8eov" << std::endl;
		return(true);
	}

	if(value == bcf_int8_missing){
		this->buffer_data_uncompressed += (S32)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int8miss" << std::endl;
		return(true);
	}

	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const S16& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added S16" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int16_vector_end){
		this->buffer_data_uncompressed += (S32)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int16eov" << std::endl;
		return(true);
	}

	if(value == bcf_int16_missing){
		this->buffer_data_uncompressed += (S32)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int16miss" << std::endl;
		return(true);
	}

	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const S32& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added S32" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int32_vector_end){
		this->buffer_data_uncompressed += (S32)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int32eov" << std::endl;
		return(true);
	}

	if(value == bcf_int32_missing){
		this->buffer_data_uncompressed += (S32)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int32miss" << std::endl;
		return(true);
	}

	this->buffer_data_uncompressed += (S32)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const U64& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = false;
	}

	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_64B, false)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added U64" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed += (U64)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const S64& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = true;
	}


	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_64B, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added S64" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed += (U64)value;
	++this->header.n_additions;
	//++this->n_entries;
	return(true);
}

bool DataContainer::Add(const float& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_FLOAT);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_FLOAT, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added FLOAT" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed += (float)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const double& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_DOUBLE);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_DOUBLE, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added DOUBLE" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed += (double)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::AddCharacter(const char& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed += (char)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::AddCharacter(const char* const string, const U32 l_string){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
		//std::cerr << "triggering: string" << std::endl;
	}

	// Make checks
	if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->buffer_data_uncompressed.Add(string, l_string);
	this->header.n_additions += l_string;
	return(true);
}

}
}
