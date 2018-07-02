#ifndef CORE_BASE_STREAMCONTAINERHEADER_H_
#define CORE_BASE_STREAMCONTAINERHEADER_H_

#include <iostream>
#include <fstream>
#include <cstring>

#include "data_container_header_controller.h"
#include "data_container_header_object.h"
#include "io/basic_buffer.h"
#include "support/enums.h"
#include "support/type_definitions.h"

namespace tachyon{
namespace containers{

struct DataContainerHeader{
private:
	typedef DataContainerHeader       self_type;
	typedef DataContainerHeaderObject header_type;
	typedef io::BasicBuffer           buffer_type;

public:
	DataContainerHeader() : identifier(0), n_entries(0), n_additions(0), n_strides(0){}
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

	// Comparators
	const bool operator==(const self_type& other) const{
		if(this->identifier    != other.identifier)    return false;
		if(this->n_entries     != other.n_entries)     return false;
		if(this->n_additions   != other.n_additions)   return false;
		if(this->n_strides     != other.n_strides)     return false;
		if(this->data_header   != other.data_header)   return false;
		if(this->stride_header != other.stride_header) return false;
		return true;
	}
	inline const bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	self_type& operator+=(const self_type& other){
		this->n_entries     += other.n_entries;
		this->n_additions   += other.n_additions;
		this->n_strides     += other.n_strides;
		return(*this);
	}

	// Accessors
	inline S32& getGlobalKey(void){ return(this->data_header.global_key); }
	inline const S32& getGlobalKey(void) const{ return(this->data_header.global_key); }
	inline const bool hasMixedStride(void) const{ return(this->data_header.hasMixedStride()); }

private:
	friend buffer_type& operator<<(buffer_type& buffer, const self_type& entry){
		buffer += entry.identifier;
		buffer += entry.n_entries;
		buffer += entry.n_additions;
		buffer += entry.n_strides;
		buffer << entry.data_header;

		if(entry.data_header.hasMixedStride()) buffer << entry.stride_header;

		return(buffer);
	}

	friend buffer_type& operator>>(buffer_type& buffer, self_type& entry){
		buffer >> entry.identifier;
		buffer >> entry.n_entries;
		buffer >> entry.n_additions;
		buffer >> entry.n_strides;
		buffer >> entry.data_header;

		if(entry.data_header.hasMixedStride()) buffer >> entry.stride_header;

		return(buffer);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.identifier),  sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.n_entries),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_additions), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_strides),   sizeof(U32));
		stream << entry.data_header;
		if(entry.data_header.hasMixedStride())
			stream << entry.stride_header;

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.identifier),  sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.n_entries),   sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_additions), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_strides),   sizeof(U32));
		stream >> entry.data_header;
		if(entry.data_header.hasMixedStride())
			stream >> entry.stride_header;

		return(stream);
	}

public:
	U64 identifier;
	U32 n_entries;      // number of container entries
	U32 n_additions;    // number of times an addition operation was executed
	U32 n_strides;      // number of stride elements
	header_type data_header;
	header_type stride_header;
};

}
}

#endif /* CORE_BASE_STREAMCONTAINERHEADER_H_ */