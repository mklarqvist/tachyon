#ifndef CORE_TOMAHAWKIMPORTWRITER_H_
#define CORE_TOMAHAWKIMPORTWRITER_H_

#include <cassert>
#include <fstream>

#include "../index/index.h"
#include "../support/type_definitions.h"

namespace tachyon {

class VariantImportWriterInterface {
private:
	typedef VariantImportWriterInterface  self_type;
	typedef index::Index                  sorted_index_type;

public:
	VariantImportWriterInterface();
	virtual ~VariantImportWriterInterface();

	void writeIndex(void);
	virtual bool open(const std::string output) =0;

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
