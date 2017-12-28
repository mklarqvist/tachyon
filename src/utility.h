#ifndef TACHYON_UTILITY_H_
#define TACHYON_UTILITY_H_

#include <iostream>
#include "support/helpers.h"
#include "support/MagicConstants.h"
#include <openssl/crypto.h>
#include <openssl/evp.h>
#include <zstd.h>

// Declare extern
std::string Tachyon::Constants::LITERAL_COMMAND_LINE;
std::string Tachyon::Constants::INTERPRETED_COMMAND;

void programMessage(const bool separator = true){
	std::cerr << "Program: " << Tachyon::Constants::PROGRAM_NAME << " " << VERSION << std::endl;
	std::cerr << "Libraries: " << Tachyon::Constants::PROGRAM_NAME << '-' << Tachyon::Constants::TACHYON_LIB_VERSION << "; "
              << SSLeay_version(SSLEAY_VERSION) << "; "
              << "ZSTD-" << ZSTD_versionString() << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/Tachyon" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if(separator) std::cerr << "----------" << std::endl;
}

void programHelp(void){
	std::cerr << "Usage: " << Tachyon::Constants::PROGRAM_NAME << " [--version] [--help] <commands> <argument>" << std::endl;
	std::cerr << "Commands: import, view" << std::endl;
}

void programHelpDetailed(void){
	programHelp();
	std::cerr <<
    "\n"
	"import       import VCF/BCF to YON\n"
    "view         YON->VCF/BCF conversion, YON subset and filter\n" << std::endl;
}

#endif /* TACHYON_UTILITY_H_ */
