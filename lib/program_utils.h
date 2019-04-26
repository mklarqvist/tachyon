#ifndef TACHYON_UTILITY_H_
#define TACHYON_UTILITY_H_

#include <iostream>

#include <htslib/hts.h>
#include <openssl/crypto.h>
#include <openssl/evp.h>
#include <zstd.h>

#include "utility.h"
#include "tachyon.h"

// These utility strings have been declared extern so that
// they can be used globally throughout the program.
// The LITERAL_COMMAND_LINE string is populated with the literal
// character input from the ABI to the program.
// The INTERPRETED_COMMAND string is the internal dump of the
// settings used for a subroutine. This is very useful for
// debugging and legacy purposes.
std::string tachyon::LITERAL_COMMAND_LINE;
std::string tachyon::INTERPRETED_COMMAND;

/**<
 * Print a standard program message to standard out. This message
 * describes the current git version, library version, and linked
 * library versions.
 * @param separator Boolean regulating the presence/absence of a dashed delimiter line after the program message.
 */
void programMessage(const bool separator = true) {
	// The OpenSSL version string generally returns a visually unappealing format
	// that includes placeholder X's. These are generally preceded by two single
	// blank spaces. If we find two consecutive blank spaces then present the
	// substring from [0, match) left-inclusive.
	std::string openssl_version_string = SSLeay_version(SSLEAY_VERSION);
	size_t match = openssl_version_string.find("  ");
	if (match != std::string::npos) openssl_version_string = openssl_version_string.substr(0, match);

	// General message for program version, git version, and linked library versions.
	std::cerr << "Program:   " << tachyon::TACHYON_PROGRAM_NAME << "-" << VERSION << " (Tools for querying and manipulating variant call data)" << std::endl;
	std::cerr << "Libraries: " << tachyon::TACHYON_PROGRAM_NAME << '-' << tachyon::TACHYON_LIB_VERSION
              << "; " << openssl_version_string
              << "; ZSTD-" << ZSTD_versionString()
			  << "; htslib " << std::string(hts_version()) << std::endl;
	std::cerr << "Contact: Marcus D. R. Klarqvist <mk819@cam.ac.uk>" << std::endl;
	std::cerr << "Documentation: https://github.com/mklarqvist/tachyon" << std::endl;
	std::cerr << "License: MIT" << std::endl;
	if (separator) std::cerr << "----------" << std::endl;
}

/**<
 * Extension of the programMessage function. First prints the
 * program message to standard out and then prints out the available
 * command names and their descriptions.
 */
void programHelp(void) {
	programMessage(false);
	std::cerr << "\nUsage: " << tachyon::TACHYON_PROGRAM_NAME << " [--version] [--help] <commands> <argument>"
    "\n\n"
    "Commands:\n"
	"import  import VCF/VCF.gz/BCF to YON\n"
    "view    convert YON->VCF/BCF/YON; provides subsetting and slicing functionality\n"
	"stats   calculate comprehensive per-sample statistics\n" << std::endl;
}

#endif /* TACHYON_UTILITY_H_ */
