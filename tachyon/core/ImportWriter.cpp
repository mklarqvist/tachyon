#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../support/MagicConstants.h"
#include "../support/helpers.h"
#include "ImportWriter.h"

namespace tachyon {

ImportWriter::ImportWriter() :
	n_blocks_written(0),
	n_variants_written(0)
{}

ImportWriter::~ImportWriter(){}

bool ImportWriter::open(const std::string output){
	if(output.size() == 0){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "No output file/file prefix provided!" << std::endl;
		return false;
	}

	this->filename = output;
	this->CheckOutputNames(output);
	this->stream.open(this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if(!SILENT){
		std::cerr << utility::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + constants::OUTPUT_SUFFIX << "..." << std::endl;
	}

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
	utility::HexToBytes(constants::TACHYON_FILE_EOF, &eof_data[0]);
	this->stream.write((char*)&eof_data[0], 32);
	this->stream.flush();
}

void ImportWriter::CheckOutputNames(const std::string& input){
	std::vector<std::string> paths = utility::filePathBaseExtension(input);
	this->basePath = paths[0];
	if(this->basePath.size() > 0)
		this->basePath += '/';

	if(paths[3].size() == constants::OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &constants::OUTPUT_SUFFIX[0], constants::OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}


} /* namespace Tachyon */
