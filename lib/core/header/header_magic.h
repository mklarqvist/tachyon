#ifndef CORE_HEADER_HEADER_MAGIC_H_
#define CORE_HEADER_HEADER_MAGIC_H_

#include <ostream>
#include <istream>
#include <cstring>

#include "io/basic_buffer.h"
#include "support/magic_constants.h"

namespace tachyon{
namespace core{

struct HeaderMagic{
public:
	typedef HeaderMagic self_type;

public:
	HeaderMagic();
	HeaderMagic(const self_type& other);
	~HeaderMagic() = default;

	inline const U64& getNumberSamples(void) const{ return(this->n_samples); }
	inline U64& getNumberSamples(void){ return(this->n_samples); }
	inline const U32& getNumberContigs(void) const{ return(this->n_contigs); }
	inline U32& getNumberContigs(void){ return(this->n_contigs); }
	inline const U16& getController(void) const{ return(this->controller); }
	inline U16& getController(void){ return(this->controller); }

	//inline bool validateMagic(void) const{ return(strncmp(&this->magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) == 0); }
	inline bool validate(void) const{
		//return(this->validateMagic() && this->n_contigs > 0 && (this->major_version > 0 || this->minor_version > 0));
		return(this->n_contigs > 0 && (this->major_version > 0 || this->minor_version > 0));
	}

	inline const bool operator!=(const self_type& other) const{ return(!(*this == other)); }
	inline const bool operator==(const self_type& other) const{
		//if(strncmp(&this->magic_string[0], &other.magic_string[0], tachyon::constants::FILE_HEADER_LENGTH) != 0) return false;
		if(this->n_samples != other.n_samples) return false;
		if(this->n_contigs != other.n_contigs) return false;
		return true;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& header){
		//stream.write(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_MAJOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_MINOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_PATCH),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&header.controller),      sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&header.n_samples),       sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&header.n_contigs),       sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.n_info_values),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.n_format_values), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.n_filter_values), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_literals),      sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_header),        sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_header_uncompressed), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		//stream.read(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		stream.read(reinterpret_cast<char*>(&header.major_version),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.minor_version),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.patch_version), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.controller),      sizeof(U16));
		stream.read(reinterpret_cast<char*>(&header.n_samples),       sizeof(U64));
		stream.read(reinterpret_cast<char*>(&header.n_contigs),       sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.n_info_values),   sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.n_format_values), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.n_filter_values), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_literals),      sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_header),        sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_header_uncompressed), sizeof(U32));
		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& header){
		//buffer.Add(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		buffer += header.major_version;
		buffer += header.minor_version;
		buffer += header.patch_version;
		buffer += header.controller;
		buffer += header.n_samples;
		buffer += header.n_contigs;
		buffer += header.n_info_values;
		buffer += header.n_format_values;
		buffer += header.n_filter_values;
		buffer += header.l_literals;
		buffer += header.l_header;
		buffer += header.l_header_uncompressed;
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& header){
		//buffer.read(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		buffer >> header.major_version;
		buffer >> header.minor_version;
		buffer >> header.patch_version;
		buffer >> header.controller;
		buffer >> header.n_samples;
		buffer >> header.n_contigs;
		buffer >> header.n_info_values;
		buffer >> header.n_format_values;
		buffer >> header.n_filter_values;
		buffer >> header.l_literals;
		buffer >> header.l_header;
		buffer >> header.l_header_uncompressed;
		return(buffer);
	}

public:
	//char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	S32  major_version;
	S32  minor_version;
	S32  patch_version;
	U16  controller;            // controller
	U64  n_samples;             // number of samples
	U32  n_contigs;             // number of contigs
	U32  n_info_values;         // number of unique info fields
	U32  n_format_values;       // number of unique format fields
	U32  n_filter_values;       // number of unique filter fields
	U32  l_literals;            // literals length
	U32  l_header;              // compressed length
	U32  l_header_uncompressed; // uncompressed length
};

}
}

#endif /* CORE_HEADER_HEADER_MAGIC_H_ */
