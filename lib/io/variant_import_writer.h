#ifndef IO_VARIANT_IMPORTWRITER_H_
#define IO_VARIANT_IMPORTWRITER_H_

#include <cassert>
#include <fstream>

#include "index/index.h"
#include "support/magic_constants.h"
#include "containers/data_container.h"

namespace tachyon {

class VariantWriterInterface {
public:
	typedef VariantWriterInterface     self_type;
	typedef index::Index               sorted_index_type;
	typedef containers::DataContainer  container_type;

public:
	VariantWriterInterface();
	virtual ~VariantWriterInterface();

	void WriteIndex(void);
	virtual bool open(const std::string output) =0;

	bool WriteBlockFooter(const container_type& footer){
		if(this->stream == nullptr) return false;
		const uint64_t start_footer_pos = this->stream->tellp();
		utility::SerializePrimitive(footer.header.data_header.uLength, *this->stream);
		utility::SerializePrimitive(footer.header.data_header.cLength, *this->stream);
		this->stream->write(reinterpret_cast<const char*>(&footer.header.data_header.crc[0]), MD5_DIGEST_LENGTH);
		*this->stream << footer.data;
		return(this->stream->good());
	}

	bool WriteEndOfBlock(void){
		if(this->stream == nullptr) return false;
		utility::SerializePrimitive(TACHYON_BLOCK_EOF, *this->stream);
		return(this->stream->good());
	}

public:
	uint64_t n_blocks_written;
	uint64_t n_variants_written;
	std::ostream* stream;
	sorted_index_type index;
};

class VariantWriterFile : public VariantWriterInterface{
public:
	typedef VariantWriterFile self_type;

public:
	VariantWriterFile();
	~VariantWriterFile();
	bool open(const std::string output);

private:
	void CheckOutputNames(const std::string& input);

public:
	// Stream information
	std::string filename;
	std::string basePath;
	std::string baseName;
};

class VariantWriterStream : public VariantWriterInterface{
public:
	typedef VariantWriterStream self_type;

public:
	VariantWriterStream();
	~VariantWriterStream();
	bool open(const std::string output){ return true; }

};

}

#endif /* IO_VARIANT_IMPORTWRITER_H_ */
