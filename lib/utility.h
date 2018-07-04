#ifndef TACHYON_UTILITY_H_
#define TACHYON_UTILITY_H_

#include <iostream>

#include <openssl/crypto.h>
#include <openssl/evp.h>
#include <zstd.h>

#include "support/helpers.h"
#include "support/magic_constants.h"

// Declare extern
std::string tachyon::constants::LITERAL_COMMAND_LINE;
std::string tachyon::constants::INTERPRETED_COMMAND;

void programMessage(const bool separator = true){
	std::cerr << "Program: " << tachyon::constants::PROGRAM_NAME << " " << VERSION << std::endl;
	std::cerr << "Libraries: " << tachyon::constants::PROGRAM_NAME << '-' << tachyon::constants::TACHYON_LIB_VERSION << "; "
              << SSLeay_version(SSLEAY_VERSION) << "; "
              << "ZSTD-" << ZSTD_versionString() << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/tachyon" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if(separator) std::cerr << "----------" << std::endl;
}

void programHelp(void){
	std::cerr << "Usage: " << tachyon::constants::PROGRAM_NAME << " [--version] [--help] <commands> <argument>" << std::endl;
	std::cerr << "Commands: import, view" << std::endl;
}

void programHelpDetailed(void){
	programHelp();
	std::cerr <<
    "\n"
	"import       import VCF/BCF to YON\n"
    "view         YON->VCF/BCF conversion, YON subset and filter\n"
	"check        comprehensive file integrity checks\n" << std::endl;
}

#endif /* TACHYON_UTILITY_H_ */
