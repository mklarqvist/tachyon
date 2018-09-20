#ifndef CORE_BASE_STREAMCONTAINERHEADER_H_
#define CORE_BASE_STREAMCONTAINERHEADER_H_

#include <iostream>
#include <fstream>
#include <cstring>

#include "data_container_header_controller.h"
#include "data_container_header_object.h"
#include "io/basic_buffer.h"

namespace tachyon{
namespace containers{

struct DataContainerHeader{
private:
	typedef DataContainerHeader       self_type;
	typedef DataContainerHeaderObject header_type;
	typedef io::BasicBuffer           buffer_type;

public:
	DataContainerHeader() :
		identifier(0),
		n_entries(0),
		n_additions(0),
		n_strides(0)
	{
	}

	DataContainerHeader(const self_type& other) :
		identifier(other.identifier),
		n_entries(other.n_entries),
		n_additions(other.n_additions),
		n_strides(other.n_strides),
		data_header(other.data_header),
		stride_header(other.stride_header)
	{

	}

	~DataContainerHeader(){}

	void reset(void){
		this->identifier  = 0;
		this->n_entries   = 0;
		this->n_additions = 0;
		this->n_strides   = 0;
		this->data_header.reset();
		this->stride_header.reset();
	}

	self_type& operator=(const self_type& other){
		this->identifier    = other.identifier;
		this->n_entries     = other.n_entries;
		this->n_additions   = other.n_additions;
		this->n_strides     = other.n_strides;
		this->data_header   = other.data_header;
		this->stride_header = other.stride_header;
		return(*this);
	}

	self_type& operator=(self_type&& other) noexcept{
		this->identifier    = other.identifier;
		this->n_entries     = other.n_entries;
		this->n_additions   = other.n_additions;
		this->n_strides     = other.n_strides;
		this->data_header   = std::move(other.data_header);
		this->stride_header = std::move(other.stride_header);
		return(*this);
	}

	// Comparators
	bool operator==(const self_type& other) const {
		if(this->identifier    != other.identifier)    return false;
		if(this->n_entries     != other.n_entries)     return false;
		if(this->n_additions   != other.n_additions)   return false;
		if(this->n_strides     != other.n_strides)     return false;
		if(this->data_header   != other.data_header)   return false;
		if(this->stride_header != other.stride_header) return false;
		return true;
	}
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	self_type& operator+=(const self_type& other) {
		this->n_entries     += other.n_entries;
		this->n_additions   += other.n_additions;
		this->n_strides     += other.n_strides;
		if(data_header.stride != other.data_header.stride || other.data_header.controller.mixedStride){
			//if(data_header.HasMixedStride() == false)
			//	std::cerr << "triggering mixed stride: " << data_header.stride << "!=" << other.data_header.stride << std::endl;
			this->data_header.SetMixedStride(true);
		}
		return(*this);
	}

	// Accessors
	inline int32_t& GetGlobalKey(void){ return(this->data_header.global_key); }
	inline const int32_t& GetGlobalKey(void) const{ return(this->data_header.global_key); }
	inline bool HasMixedStride(void) const{ return(this->data_header.HasMixedStride()); }

private:
	friend buffer_type& operator<<(buffer_type& buffer, const self_type& entry){
		buffer += entry.identifier;
		buffer += entry.n_entries;
		buffer += entry.n_additions;
		buffer += entry.n_strides;
		buffer << entry.data_header;

		if(entry.data_header.HasMixedStride()) buffer << entry.stride_header;

		return(buffer);
	}

	friend buffer_type& operator>>(buffer_type& buffer, self_type& entry){
		buffer >> entry.identifier;
		buffer >> entry.n_entries;
		buffer >> entry.n_additions;
		buffer >> entry.n_strides;
		buffer >> entry.data_header;

		if(entry.data_header.HasMixedStride()) buffer >> entry.stride_header;

		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.identifier),  sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_entries),   sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_additions), sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_strides),   sizeof(uint32_t));
		stream << entry.data_header;
		if(entry.data_header.HasMixedStride())
			stream << entry.stride_header;

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.identifier),  sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.n_entries),   sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_additions), sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_strides),   sizeof(uint32_t));
		stream >> entry.data_header;
		if(entry.data_header.HasMixedStride())
			stream >> entry.stride_header;

		return(stream);
	}

public:
	uint64_t identifier;
	uint32_t n_entries;      // number of container entries
	uint32_t n_additions;    // number of times an addition operation was executed
	uint32_t n_strides;      // number of stride elements
	header_type data_header;
	header_type stride_header;
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADER_H_ */
