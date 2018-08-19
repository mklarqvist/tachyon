#ifndef CORE_TOMAHAWKIMPORTWRITER_H_
#define CORE_TOMAHAWKIMPORTWRITER_H_

#include <cassert>
#include <fstream>

#include "index/index.h"
#include "support/magic_constants.h"
#include "containers/data_container.h"

namespace tachyon {

class VariantImportWriterInterface {
private:
	typedef VariantImportWriterInterface  self_type;
	typedef index::Index                  sorted_index_type;
	typedef containers::DataContainer     container_type;

public:
	VariantImportWriterInterface();
	virtual ~VariantImportWriterInterface();

	void writeIndex(void);
	virtual bool open(const std::string output) =0;

	bool WriteBlockFooter(const container_type& footer){
		if(this->stream == nullptr) return false;
		const uint64_t start_footer_pos = this->stream->tellp();
		utility::SerializePrimitive(footer.header.data_header.uLength, *this->stream);
		utility::SerializePrimitive(footer.header.data_header.cLength, *this->stream);
		this->stream->write(reinterpret_cast<const char*>(&footer.header.data_header.crc[0]), MD5_DIGEST_LENGTH);
		*this->stream << footer.buffer_data;
		return(this->stream->good());
	}

	bool WriteEndOfBlock(void){
		if(this->stream == nullptr) return false;
		utility::SerializePrimitive(constants::TACHYON_BLOCK_EOF, *this->stream);
		return(this->stream->good());
	}

public:
	uint64_t n_blocks_written;
	uint64_t n_variants_written;
	std::ostream* stream;
	sorted_index_type index;
};

class VariantImportWriterFile : public VariantImportWriterInterface{
private:
	typedef VariantImportWriterFile self_type;

public:
	VariantImportWriterFile();
	~VariantImportWriterFile();
	bool open(const std::string output);

private:
	void checkOutputNames(const std::string& input);

public:
	// Stream information
	std::string   filename;
	std::string   basePath;
	std::string   baseName;
};

class VariantImportWriterStream : public VariantImportWriterInterface{
private:
	typedef VariantImportWriterStream self_type;

public:
	VariantImportWriterStream();
	~VariantImportWriterStream();
	bool open(const std::string output){ return true; }

};

} /* namespace Tomahawk */

#endif /* CORE_TOMAHAWKIMPORTWRITER_H_ */
