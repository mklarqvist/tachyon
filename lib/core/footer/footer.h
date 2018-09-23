#ifndef CORE_FOOTER_FOOTER_H_
#define CORE_FOOTER_FOOTER_H_

#include "support/magic_constants.h"
#include "support/helpers.h"

namespace tachyon{
namespace core{

#define YON_FOOTER_LENGTH ((TACHYON_FILE_EOF_LENGTH) + sizeof(uint64_t)*3 + sizeof(uint16_t))

struct Footer {
public:
	typedef Footer self_type;

public:
	Footer() :
		offset_end_of_data(0),
		n_blocks(0),
		n_variants(0),
		controller(0)
	{
		utility::HexToBytes(TACHYON_FILE_EOF, &this->EOF_marker[0]);
	}

	Footer(const char* const data) :
		offset_end_of_data(*reinterpret_cast<const uint64_t* const>(data)),
		n_blocks(*reinterpret_cast<const uint64_t* const>(&data[sizeof(uint64_t)])),
		n_variants(*reinterpret_cast<const uint64_t* const>(&data[sizeof(uint64_t)*2])),
		controller(*reinterpret_cast<const uint16_t* const>(&data[sizeof(uint64_t)*3]))
	{
		memcpy(&this->EOF_marker[0], &data[sizeof(uint64_t)*3+sizeof(uint16_t)], TACHYON_FILE_EOF_LENGTH);
	}

	Footer(const self_type& other) :
		offset_end_of_data(other.offset_end_of_data),
		n_blocks(other.n_blocks),
		n_variants(other.n_variants),
		controller(other.controller)
	{
		memcpy(&this->EOF_marker[0], &other.EOF_marker[0], TACHYON_FILE_EOF_LENGTH);
	}

	~Footer() = default;

	inline const uint64_t& GetEODOffset(void) const{ return(this->offset_end_of_data); }
	inline uint64_t& GetEODOffset(void){ return(this->offset_end_of_data); }
	inline const uint64_t& GetNumberBlocks(void) const{ return(this->n_blocks); }
	inline uint64_t& GetNumberBlocks(void){ return(this->n_blocks); }
	inline const uint64_t& GetNumberVariants(void) const{ return(this->n_variants); }
	inline uint64_t& GetNumberVariants(void){ return(this->n_variants); }
	inline const uint16_t& GetController(void) const{ return(this->controller); }
	inline uint16_t& GetController(void){ return(this->controller); }

	inline bool Validate(void) const{
		if(this->offset_end_of_data == 0) return false;
		if(this->n_blocks  == 0)          return false;
		if(this->n_variants == 0)         return false;

		// Check EOF marker
		uint8_t reference[TACHYON_FILE_EOF_LENGTH];
		utility::HexToBytes(TACHYON_FILE_EOF, &reference[0]);

		if(strncmp(reinterpret_cast<const char* const>(&this->EOF_marker[0]), reinterpret_cast<const char* const>(&reference[0]), TACHYON_FILE_EOF_LENGTH) != 0) return false;
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& footer){
		stream.write(reinterpret_cast<const char*>(&footer.offset_end_of_data), sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&footer.n_blocks),           sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&footer.n_variants),         sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&footer.controller),         sizeof(uint16_t));
		stream.write(reinterpret_cast<const char*>(&footer.EOF_marker[0]), TACHYON_FILE_EOF_LENGTH);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& footer){
		stream.read(reinterpret_cast<char*>(&footer.offset_end_of_data), sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&footer.n_blocks),           sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&footer.n_variants),         sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&footer.controller),         sizeof(uint16_t));
		stream.read(reinterpret_cast<char*>(&footer.EOF_marker[0]), TACHYON_FILE_EOF_LENGTH);
		return(stream);
	}

public:
	uint64_t  offset_end_of_data;
	uint64_t  n_blocks;
	uint64_t  n_variants;
	uint16_t  controller;
    uint8_t   EOF_marker[TACHYON_FILE_EOF_LENGTH];
};

}
}

#endif /* CORE_FOOTER_FOOTER_H_ */
