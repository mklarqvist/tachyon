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
		const U64 start_footer_pos = this->stream->tellp();
		const U32 footer_uLength   = footer.header.data_header.uLength;
		const U32 footer_cLength   = footer.header.data_header.cLength;
		const U32 footer_crc       = footer.header.data_header.crc;
		this->stream->write(reinterpret_cast<const char*>(&footer_uLength), sizeof(U32));
		this->stream->write(reinterpret_cast<const char*>(&footer_cLength), sizeof(U32));
		this->stream->write(reinterpret_cast<const char*>(&footer_crc),     sizeof(U32));
		*this->stream << footer.buffer_data;
		return(this->stream->good());
	}

	bool WriteEndOfBlock(void){
		if(this->stream == nullptr) return false;
		this->stream->write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));
		return(this->stream->good());
	}

public:
	U64 n_blocks_written;
	U64 n_variants_written;
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
