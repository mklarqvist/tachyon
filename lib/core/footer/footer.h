#ifndef CORE_FOOTER_FOOTER_H_
#define CORE_FOOTER_FOOTER_H_


#include "support/type_definitions.h"
#include "support/magic_constants.h"
#include "support/helpers.h"

namespace tachyon{
namespace core{

#define YON_FOOTER_LENGTH ((constants::TACHYON_FILE_EOF_LENGTH) + sizeof(U64)*3 + sizeof(U16))

struct Footer{
public:
	typedef Footer self_type;

public:
	Footer() :
		offset_end_of_data(0),
		n_blocks(0),
		n_variants(0),
		controller(0)
	{
		utility::HexToBytes(constants::TACHYON_FILE_EOF, &this->EOF_marker[0]);
	}

	Footer(const char* const data) :
		offset_end_of_data(*reinterpret_cast<const U64* const>(data)),
		n_blocks(*reinterpret_cast<const U64* const>(&data[sizeof(U64)])),
		n_variants(*reinterpret_cast<const U64* const>(&data[sizeof(U64)*2])),
		controller(*reinterpret_cast<const U16* const>(&data[sizeof(U64)*3]))
	{
		memcpy(&this->EOF_marker[0], &data[sizeof(U64)*3+sizeof(U16)], constants::TACHYON_FILE_EOF_LENGTH);
	}

	Footer(const self_type& other) :
		offset_end_of_data(other.offset_end_of_data),
		n_blocks(other.n_blocks),
		n_variants(other.n_variants),
		controller(other.controller)
	{
		memcpy(&this->EOF_marker[0], &other.EOF_marker[0], constants::TACHYON_FILE_EOF_LENGTH);
	}

	~Footer() = default;

	inline const U64& getEODOffset(void) const{ return(this->offset_end_of_data); }
	inline U64& getEODOffset(void){ return(this->offset_end_of_data); }
	inline const U64& getNumberBlocks(void) const{ return(this->n_blocks); }
	inline U64& getNumberBlocks(void){ return(this->n_blocks); }
	inline const U64& getNumberVariants(void) const{ return(this->n_variants); }
	inline U64& getNumberVariants(void){ return(this->n_variants); }
	inline const U16& getController(void) const{ return(this->controller); }
	inline U16& getController(void){ return(this->controller); }

	inline const bool validate(void) const{
		if(this->offset_end_of_data == 0) return false;
		if(this->n_blocks  == 0)          return false;
		if(this->n_variants == 0)         return false;

		// Check EOF marker
		BYTE reference[constants::TACHYON_FILE_EOF_LENGTH];
		utility::HexToBytes(constants::TACHYON_FILE_EOF, &reference[0]);

		if(strncmp(reinterpret_cast<const char* const>(&this->EOF_marker[0]), reinterpret_cast<const char* const>(&reference[0]), constants::TACHYON_FILE_EOF_LENGTH) != 0) return false;
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& footer){
		stream.write(reinterpret_cast<const char*>(&footer.offset_end_of_data), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&footer.n_blocks),           sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&footer.n_variants),         sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&footer.controller),         sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&footer.EOF_marker[0]), constants::TACHYON_FILE_EOF_LENGTH);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& footer){
		stream.read(reinterpret_cast<char*>(&footer.offset_end_of_data), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&footer.n_blocks),           sizeof(U64));
		stream.read(reinterpret_cast<char*>(&footer.n_variants),         sizeof(U64));
		stream.read(reinterpret_cast<char*>(&footer.controller),         sizeof(U16));
		stream.read(reinterpret_cast<char*>(&footer.EOF_marker[0]), constants::TACHYON_FILE_EOF_LENGTH);
		return(stream);
	}

public:
	U64  offset_end_of_data;
	U64  n_blocks;
	U64  n_variants;
	U16  controller;
    BYTE EOF_marker[constants::TACHYON_FILE_EOF_LENGTH];
};

}
}

#endif /* CORE_FOOTER_FOOTER_H_ */
