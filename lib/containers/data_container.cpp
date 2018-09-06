#include <cassert>

#include "data_container.h"
#include "io/basic_buffer.h"
#include "support/helpers.h"
#include "algorithm/digest/variant_digest_manager.h"

namespace tachyon{
namespace containers{

DataContainer::DataContainer()
{}

DataContainer::DataContainer(const uint32_t start_size) :
	data(start_size),
	strides(start_size),
	data_uncompressed(start_size),
	strides_uncompressed(start_size)
{}

DataContainer::~DataContainer(){ }

DataContainer::DataContainer(self_type&& other) noexcept :
	header(std::move(other.header)),
	data(std::move(other.data)),
	strides(std::move(other.strides)),
	data_uncompressed(std::move(other.data_uncompressed)),
	strides_uncompressed(std::move(other.strides_uncompressed))
{

}

DataContainer::DataContainer(const self_type& other) :
	header(other.header),
	data(other.data),
	strides(other.strides),
	data_uncompressed(other.data_uncompressed),
	strides_uncompressed(other.strides_uncompressed)
{}

DataContainer& DataContainer::operator=(const self_type& other){
	this->data = other.data;
	this->data_uncompressed = other.data_uncompressed;
	this->strides = other.strides;
	this->strides_uncompressed = other.strides_uncompressed;
	this->header = other.header;
	return(*this);
}

DataContainer& DataContainer::operator=(self_type&& other) noexcept{
	this->data = std::move(other.data);
	this->data_uncompressed = std::move(other.data_uncompressed);
	this->strides = std::move(other.strides);
	this->strides_uncompressed = std::move(other.strides_uncompressed);
	this->header = std::move(other.header);
	return(*this);
}

void DataContainer::reset(void){
	this->data.reset();
	this->data_uncompressed.reset();
	this->strides.reset();
	this->strides_uncompressed.reset();
	this->header.reset();
}

void DataContainer::resize(const uint32_t size){
	this->data.resize(size);
	this->data_uncompressed.resize(size);
	this->strides.resize(size);
	this->strides_uncompressed.resize(size);
}

void DataContainer::GenerateMd5(void){
	algorithm::VariantDigestManager::GenerateMd5(this->data_uncompressed.data(), this->data_uncompressed.size(), &this->header.data_header.crc[0]);
	algorithm::VariantDigestManager::GenerateMd5(this->strides_uncompressed.data(), this->strides_uncompressed.size(), &this->header.stride_header.crc[0]);
}

bool DataContainer::CheckMd5(int target){
	if(target == 0){
		if(this->data_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->data_uncompressed.data(), this->data_uncompressed.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if(target == 1){
		if(this->strides_uncompressed.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->strides_uncompressed.data(), this->strides_uncompressed.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));

	} else if(target == 3){
		if(this->data.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->data.data(), this->data.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if(target == 4){
		if(this->strides.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->strides.data(), this->strides.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));
	}
	return true;
}

bool DataContainer::CheckUniformity(void){
	if(data_uncompressed.size() == 0)
		return false;

	if(this->header.n_entries == 0)
		return false;

	// We know the stride cannot be uniform if
	// the stride size is uneven
	const int16_t& stride_size = this->header.data_header.stride;
	if(stride_size == -1)
		return false;

	uint32_t stride_update = stride_size;

	uint8_t word_width = sizeof(char);
	switch(this->header.data_header.controller.type){
	case YON_TYPE_DOUBLE: stride_update *= sizeof(double);   word_width = sizeof(double);  break;
	case YON_TYPE_FLOAT:  stride_update *= sizeof(float);    word_width = sizeof(float);   break;
	case YON_TYPE_8B:     stride_update *= sizeof(uint8_t);  word_width = sizeof(uint8_t); break;
	case YON_TYPE_16B:    stride_update *= sizeof(uint16_t); word_width = sizeof(uint16_t);break;
	case YON_TYPE_32B:    stride_update *= sizeof(int32_t);  word_width = sizeof(int32_t); break;
	case YON_TYPE_64B:    stride_update *= sizeof(uint64_t); word_width = sizeof(uint64_t);break;
	case YON_TYPE_CHAR:   stride_update *= sizeof(char);     word_width = sizeof(char);    break;
	default: return false; break;
	}

	const uint64_t first_hash = XXH64(this->data_uncompressed.data(), stride_update, 2147483647);

	uint64_t cumulative_position = stride_update;
	for(uint32_t i = 1; i < this->header.n_entries; ++i){
		if(XXH64(&this->data_uncompressed.buffer_[cumulative_position], stride_update, 2147483647) != first_hash){
			return(false);
		}
		cumulative_position += stride_update;
	}
	assert(cumulative_position == this->data_uncompressed.size());

	this->header.n_entries   = 1;
	this->header.n_strides   = 0;
	this->data_uncompressed.n_chars_                = stride_size * word_width;
	this->header.data_header.uLength                = stride_size * word_width;
	this->header.data_header.cLength                = stride_size * word_width;
	this->header.data_header.controller.uniform     = true;
	this->header.data_header.controller.mixedStride = false;
	this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
	return(true);
}

void DataContainer::ReformatInteger(){
	if(data_uncompressed.size() == 0)
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
	assert(this->data_uncompressed.size() % sizeof(int32_t) == 0);
	assert(this->header.n_additions * sizeof(int32_t) == this->data_uncompressed.size());
	const int32_t* const dat  = reinterpret_cast<const int32_t* const>(this->data_uncompressed.data());
	int32_t min_value = dat[0];
	int32_t max_value = dat[0];
	bool has_special = false;

	// Iterate over available data and search for either missingness
	// or sentinel node values. If a match is found trigger the
	// bool flag to signal the use of the recoding procedure.
	for(uint32_t j = 0; j < this->header.n_additions; ++j){
		if(dat[j] == bcf_int32_missing || dat[j] == bcf_int32_vector_end){
			has_special = true;
			continue;
		}
		if(dat[j] < min_value) min_value = dat[j];
		if(dat[j] > max_value) max_value = dat[j];
	}

	// If we have missing values then we have to use signed
	// primitives to accommodate this fact.
	uint8_t byte_width = 0;
	if(min_value < 0 || has_special == true){
		byte_width = ceil((ceil(log2(abs(min_value) + 1 + 2)) + 1) / 8);  // One bit is used for sign, 2 values for missing and end-of-vector
		const uint8_t byte_width_max = ceil((ceil(log2(abs(max_value) + 1 + 2)) + 1) / 8);
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
	this->data.reset();
	this->data.resize(this->data_uncompressed.size() + 65536);

	// Is non-negative
	// Also cannot have missing values
	if(min_value >= 0 && has_special == false){
		this->header.data_header.controller.signedness = 0;

		if(byte_width == 1){
			this->header.data_header.controller.type = YON_TYPE_8B;

			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				this->data += (uint8_t)dat[j];
			}

		} else if(byte_width == 2){
			this->header.data_header.controller.type = YON_TYPE_16B;

			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				this->data += (uint16_t)dat[j];
			}

		} else if(byte_width == 4){
			this->header.data_header.controller.type = YON_TYPE_32B;

			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				this->data += (uint32_t)dat[j];
			}

		} else if(byte_width == 8){
			this->header.data_header.controller.type = YON_TYPE_64B;

			for(uint32_t j = 0; j < this->header.n_additions; ++j)
				this->data += (uint64_t)dat[j];

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

			const int8_t missing = INT8_MIN;
			const int8_t eov     = INT8_MIN + 1;
			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int8_t)dat[j];
			}

		} else if(byte_width == 2){
			this->header.data_header.controller.type = YON_TYPE_16B;

			const int16_t missing = INT16_MIN;
			const int16_t eov     = INT16_MIN + 1;
			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int16_t)dat[j];
			}

		} else if(byte_width == 4){
			this->header.data_header.controller.type = YON_TYPE_32B;

			const int32_t missing = INT32_MIN;
			const int32_t eov     = INT32_MIN + 1;
			for(uint32_t j = 0; j < this->header.n_additions; ++j){
				if(dat[j] == bcf_int32_missing)         this->data += missing;
				else if(dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int32_t)dat[j];
			}

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}
	assert(this->data.size() % byte_width == 0);
	assert(this->header.n_additions * byte_width == this->data.size());

	memcpy(this->data_uncompressed.data(), this->data.data(), this->data.size());
	this->data_uncompressed.n_chars_ = this->data.size();
	this->header.data_header.uLength        = this->data_uncompressed.size();
	this->data.reset();
}

void DataContainer::ReformatStride(){
	if(this->strides_uncompressed.size() == 0)
		return;

	if(this->header.data_header.HasMixedStride() == false)
		return;

	// Recode integer types
	if(!(this->header.stride_header.controller.type == YON_TYPE_32B &&
	   this->header.stride_header.controller.signedness == 0)){
		std::cerr << utility::timestamp("ERROR") << "Illegal to have non-uint32_t values at this point: " << this->header.stride_header.controller.type << ":" << this->header.stride_header.controller.signedness << std::endl;
		exit(1);
	}

	// At this point all integers are uint32_t
	const uint32_t* const dat = reinterpret_cast<const uint32_t* const>(this->strides_uncompressed.data());
	assert(this->strides_uncompressed.size() % sizeof(uint32_t) == 0);
	assert(this->strides_uncompressed.size() / sizeof(uint32_t) == this->header.n_strides);

	uint32_t max = 1;
	for(uint32_t j = 0; j < this->header.n_strides; ++j){
		if(dat[j] > max)
			max = dat[j];
	}

	uint8_t byte_width = ceil(log2(max + 1) / 8);
	if(byte_width == 0) byte_width = 1;
	if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
	if(byte_width > 4)  byte_width = 8;

	// This cannot ever be uniform
	this->strides.reset();
	this->strides.resize(this->strides_uncompressed.size() + 65536);

	if(byte_width == 1){
		this->header.stride_header.controller.type = YON_TYPE_8B;

		for(uint32_t j = 0; j < this->header.n_strides; ++j){
			//assert((uint8_t)dat[j] == dat[j]);
			this->strides += (uint8_t)dat[j];
		}

	} else if(byte_width == 2){
		this->header.stride_header.controller.type = YON_TYPE_16B;

		for(uint32_t j = 0; j < this->header.n_strides; ++j){
			//assert((uint16_t)dat[j] == dat[j]);
			this->strides += (uint16_t)dat[j];
		}

	} else if(byte_width == 4){
		this->header.stride_header.controller.type = YON_TYPE_32B;

		for(uint32_t j = 0; j < this->header.n_strides; ++j){
			//assert((uint32_t)dat[j] == dat[j]);
			this->strides += (uint32_t)dat[j];
		}

	} else if(byte_width == 8){
		this->header.stride_header.controller.type = YON_TYPE_64B;

		for(uint32_t j = 0; j < this->header.n_strides; ++j){
			//assert((uint64_t)dat[j] == dat[j]);
			this->strides += (uint64_t)dat[j];
		}

	} else {
		std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
		exit(1);
	}
	memcpy(this->strides_uncompressed.data(), this->strides.data(), this->strides.size());
	this->strides_uncompressed.n_chars_ = this->strides.size();
	this->header.stride_header.uLength         = this->strides_uncompressed.size();
	this->strides.reset();
}

uint32_t DataContainer::GetObjectSize(void) const{
	// In case data is encrypted
	if(this->header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
		return(this->data.size());

	uint32_t total_size = this->data.size();
	if(this->header.data_header.HasMixedStride())
		total_size += this->strides.size();

	return(total_size);
}

uint64_t DataContainer::GetObjectSizeUncompressed(void) const{
	uint64_t total_size = this->data_uncompressed.size();
	if(this->header.data_header.HasMixedStride())
		total_size += this->strides_uncompressed.size();

	return(total_size);
}

void DataContainer::UpdateContainer(bool reformat_data, bool reformat_stride){
	// If the data container Has entries in it but has
	// no actual data then it is a BOOLEAN
	if(this->header.n_entries && this->data_uncompressed.size() == 0){
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

	if(this->data_uncompressed.size() == 0)
		return;

	// Check if stream is uniform in content
	if(this->header.data_header.controller.type != YON_TYPE_STRUCT){
		this->CheckUniformity();
		// Reformat stream to use as small word size as possible
		if(reformat_data) this->ReformatInteger();
	}

	// Set uncompressed length
	this->header.data_header.uLength = this->data_uncompressed.size();

	// If we have mixed striding
	if(this->header.data_header.HasMixedStride()){
		// Reformat stream to use as small word size as possible
		if(reformat_stride) this->ReformatStride();
		this->header.stride_header.uLength = this->strides_uncompressed.size();
	}
}

void DataContainer::AddStride(const uint32_t value){
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
	this->strides_uncompressed += (uint32_t)value;
	++this->header.n_strides;
}

bool DataContainer::Add(const uint8_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint8_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const uint16_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint16_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const uint32_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()){
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint32_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const int8_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}

	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int8_t" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int8_vector_end){
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int8eov" << std::endl;
		return(true);
	}

	if(value == bcf_int8_missing){
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int8miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const int16_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int16_t" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int16_vector_end){
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int16eov" << std::endl;
		return(true);
	}

	if(value == bcf_int16_missing){
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int16miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const int32_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if(!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int32_t" << std::endl;
		exit(1);
		return false;
	}

	if(value == bcf_int32_vector_end){
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int32eov" << std::endl;
		return(true);
	}

	if(value == bcf_int32_missing){
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int32miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const uint64_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = false;
	}

	// Make checks
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_64B, false)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint64_t" << std::endl;
		return false;
	}

	this->data_uncompressed += (uint64_t)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const int64_t& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = true;
	}


	// Make checks
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_64B, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int64_t" << std::endl;
		return false;
	}

	this->data_uncompressed += (uint64_t)value;
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
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_FLOAT, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added FLOAT" << std::endl;
		return false;
	}

	this->data_uncompressed += (float)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::Add(const double& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_DOUBLE);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_DOUBLE, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added DOUBLE" << std::endl;
		return false;
	}

	this->data_uncompressed += (double)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::AddCharacter(const char& value){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->data_uncompressed += (char)value;
	++this->header.n_additions;
	return(true);
}

bool DataContainer::AddCharacter(const char* const string, const uint32_t l_string){
	if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
		//std::cerr << "triggering: string" << std::endl;
	}

	// Make checks
	if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->data_uncompressed.Add(string, l_string);
	this->header.n_additions += l_string;
	return(true);
}

std::ostream& operator<<(std::ostream& stream, const DataContainer& entry){
	stream << entry.data;
	if(entry.header.data_header.HasMixedStride())
		stream << entry.strides;

	return(stream);
}

std::istream& operator>>(std::istream& stream, DataContainer& entry){
	if(entry.header.data_header.controller.encryption == YON_ENCRYPTION_NONE){
		entry.data.resize(entry.header.data_header.cLength);
		stream.read(entry.data.data(), entry.header.data_header.cLength);
		entry.data.n_chars_ = entry.header.data_header.cLength;

		if(entry.header.data_header.HasMixedStride()){
			entry.strides.resize(entry.header.stride_header.cLength);
			stream.read(entry.strides.data(), entry.header.stride_header.cLength);
			entry.strides.n_chars_ = entry.header.stride_header.cLength;
		}
	} else { // Data is encrypted
		entry.data.resize(entry.header.data_header.eLength);
		stream.read(entry.data.data(), entry.header.data_header.eLength);
		entry.data.n_chars_ = entry.header.data_header.eLength;
	}
	return(stream);
}

}
}
