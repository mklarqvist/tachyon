#ifndef CORE_HEADER_READ_HEADER_MAGIC_H_
#define CORE_HEADER_READ_HEADER_MAGIC_H_

#include <iostream>

namespace tachyon{
namespace core{

struct ReadHeaderMagic{
public:
	typedef ReadHeaderMagic self_type;

public:
	ReadHeaderMagic();
	ReadHeaderMagic(const self_type& other);
	~ReadHeaderMagic() = default;

	inline const U16& getController(void) const{ return(this->controller); }
	inline U16& getController(void){ return(this->controller); }

	inline bool validateMagic(void) const{ return(strncmp(&this->magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) == 0); }
	inline bool validate(void) const{
		return(this->validateMagic() && (this->major_version > 0 || this->minor_version > 0));
	}


private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& header){
		stream.write(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_MAJOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_MINOR),   sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&tachyon::constants::TACHYON_VERSION_RELEASE), sizeof(S32));
		stream.write(reinterpret_cast<const char*>(&header.controller),      sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&header.l_literals),      sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_header),        sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&header.l_header_uncompressed), sizeof(U32));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.magic_string, tachyon::constants::FILE_HEADER_LENGTH);
		stream.read(reinterpret_cast<char*>(&header.major_version),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.minor_version),   sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.release_version), sizeof(S32));
		stream.read(reinterpret_cast<char*>(&header.controller),      sizeof(U16));
		stream.read(reinterpret_cast<char*>(&header.l_literals),      sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_header),        sizeof(U32));
		stream.read(reinterpret_cast<char*>(&header.l_header_uncompressed), sizeof(U32));
		return(stream);
	}

public:
	char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	S32  major_version;
	S32  minor_version;
	S32  release_version;
	U16  controller;            // controller
	U32  l_literals;            // literals length
	U32  l_header;              // compressed length
	U32  l_header_uncompressed; // uncompressed length
};

}
}



#endif /* CORE_HEADER_READ_HEADER_MAGIC_H_ */
