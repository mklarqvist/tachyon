#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "variant_writer.h"
#include "utility.h"
#include "tachyon.h"

namespace tachyon {

VariantWriterInterface::VariantWriterInterface() :
	n_blocks_written(0),
	n_variants_written(0),
	stream(nullptr)
{}

VariantWriterInterface::~VariantWriterInterface()
{
}

void VariantWriterInterface::WriteIndex(void){
	*this->stream << this->index;
	this->stream->flush();
}

bool VariantWriterFile::open(const std::string output){
	if(output.size() == 0){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "No output file/file prefix provided!" << std::endl;
		return false;
	}
	std::ofstream* ostream = reinterpret_cast<std::ofstream*>(this->stream);

	this->filename = output;
	this->CheckOutputNames(output);
	ostream->open(this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->stream->good()){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << utility::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

void VariantWriterFile::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = utility::FilePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == TACHYON_OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &TACHYON_OUTPUT_SUFFIX[0], TACHYON_OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}

VariantWriterFile::VariantWriterFile(){ this->stream = new std::ofstream; }

VariantWriterFile::~VariantWriterFile(){
	this->stream->flush();
	delete this->stream;
}

VariantWriterStream::VariantWriterStream(){ this->stream = &std::cout; }

VariantWriterStream::~VariantWriterStream(){
	this->stream->flush();;
}

} /* namespace Tachyon */
