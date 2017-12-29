#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../support/MagicConstants.h"
#include "../support/helpers.h"
#include "ImportWriter.h"

namespace Tachyon {

ImportWriter::ImportWriter()
{}

ImportWriter::~ImportWriter(){}

bool ImportWriter::Open(const std::string output){
	this->filename = output;
	this->CheckOutputNames(output);
	this->stream.open(this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + Constants::OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

bool ImportWriter::WriteHeader(void){
	if(!this->stream.good()){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Stream is bad!" << std::endl;
		return false;
	}

	this->stream.write(&Constants::FILE_HEADER[0], Constants::FILE_HEADER.size());
	this->stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_VERSION_MAJOR),  sizeof(U16));
	this->stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_VERSION_MINOR),  sizeof(U16));
	this->stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_VERSION_RELEASE),sizeof(U16));
	return true;
}

void ImportWriter::WriteFinal(void){
	this->index.buildSuperIndex();
	this->stream << this->index;
}

void ImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = Helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == Constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &Constants::OUTPUT_SUFFIX[0], Constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tachyon */
