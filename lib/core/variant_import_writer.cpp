#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "variant_import_writer.h"
#include "support/helpers.h"
#include "support/magic_constants.h"

namespace tachyon {

VariantImportWriterInterface::VariantImportWriterInterface() :
	n_blocks_written(0),
	n_variants_written(0),
	stream(nullptr)
{}

VariantImportWriterInterface::~VariantImportWriterInterface()
{
}

void VariantImportWriterInterface::WriteIndex(void){
	*this->stream << this->index;
	this->stream->flush();
}

bool VariantImportWriterFile::open(const std::string output){
	if(output.size() == 0){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "No output file/file prefix provided!" << std::endl;
		return false;
	}
	std::ofstream* ostream = reinterpret_cast<std::ofstream*>(this->stream);

	this->filename = output;
	this->checkOutputNames(output);
	ostream->open(this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->stream->good()){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << utility::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

void VariantImportWriterFile::checkOutputNames(const std::string& input){
	std::vector<std::string> paths = utility::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &constants::OUTPUT_SUFFIX[0], constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}

VariantImportWriterFile::VariantImportWriterFile(){ this->stream = new std::ofstream; }
VariantImportWriterFile::~VariantImportWriterFile(){
	this->stream->flush();
	delete this->stream;
}
VariantImportWriterStream::VariantImportWriterStream(){ this->stream = &std::cout; }
VariantImportWriterStream::~VariantImportWriterStream(){
	this->stream->flush();;
}

} /* namespace Tachyon */
