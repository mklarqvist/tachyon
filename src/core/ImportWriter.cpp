#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../support/MagicConstants.h"
#include "../support/helpers.h"
#include "ImportWriter.h"

namespace tachyon {

ImportWriter::ImportWriter()
{}

ImportWriter::~ImportWriter(){}

bool ImportWriter::Open(const std::string output){
	this->filename = output;
	this->CheckOutputNames(output);
	this->stream.open(this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->stream.good()){
		std::cerr << helpers::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << helpers::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

bool ImportWriter::WriteHeader(void){
	if(!this->stream.good()){
		std::cerr << helpers::timestamp("ERROR", "WRITER") << "Stream is bad!" << std::endl;
		return false;
	}

	this->stream.write(&constants::FILE_HEADER[0], constants::FILE_HEADER.size());
	this->stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MAJOR),  sizeof(U16));
	this->stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_MINOR),  sizeof(U16));
	this->stream.write(reinterpret_cast<const char*>(&constants::TACHYON_VERSION_RELEASE),sizeof(U16));
	return true;
}

void ImportWriter::WriteIndex(void){
	this->index.buildSuperIndex();
	this->stream << this->index;
	this->stream.flush();
}

void ImportWriter::WriteFinal(const U64& data_ends){
	this->stream.write(reinterpret_cast<const char* const>(&data_ends), sizeof(U64));

	// Write EOF
	BYTE eof_data[32];
	helpers::HexToBytes(constants::TACHYON_FILE_EOF, &eof_data[0]);
	this->stream.write((char*)&eof_data[0], 32);
	this->stream.flush();
}

void ImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = helpers::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &constants::OUTPUT_SUFFIX[0], constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tachyon */
