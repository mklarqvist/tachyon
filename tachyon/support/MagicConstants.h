#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>
#include <regex>

#include "../support/type_definitions.h"

extern int SILENT;

namespace tachyon{
namespace constants{

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

/*------   Version   ------*/
const S32 TACHYON_VERSION_MAJOR    = 0;
const S32 TACHYON_VERSION_MINOR    = 1;
const S32 TACHYON_VERSION_RELEASE  = 0;
const S32 TACHYON_VERSION_NUMBER   = (TACHYON_VERSION_MAJOR *100*100 + TACHYON_VERSION_MINOR *100 + TACHYON_VERSION_RELEASE);
const std::string TACHYON_LIB_VERSION = std::to_string(TACHYON_VERSION_MAJOR) + '.' + std::to_string(TACHYON_VERSION_MINOR) + '.' + std::to_string(TACHYON_VERSION_RELEASE);

/*------   Basics   ------*/
const std::string PROGRAM_NAME  = "tachyon";
const std::string OUTPUT_SUFFIX = "yon";
const std::string FILE_HEADER   = "TACHYON\1";
const U32 FILE_HEADER_LENGTH    = 8;

/*------   Regular expression patterns  ------*/
const std::regex YON_REGEX_CONTIG_ONLY     = std::regex("^[A-Za-z0-9\\-_]+$");
const std::regex YON_REGEX_CONTIG_POSITION = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([eE]{1}[0-9]{1})?$");
const std::regex YON_REGEX_CONTIG_RANGE    = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([eE]{1}[0-9]{1})?\\-[0-9]+([eE]{1}[0-9]{1})?$");
const std::regex YON_REGEX_PACKED_ALLELES  = std::regex("^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$");

/*------   Hot meta biallelic packing   ------*/
// Encoding for bases
const char* const REF_ALT_LOOKUP = "ATGC.XN";
const BYTE REF_ALT_A = 0;
const BYTE REF_ALT_T = 1;
const BYTE REF_ALT_G = 2;
const BYTE REF_ALT_C = 3;
const BYTE REF_ALT_MISSING = 4;
const BYTE REF_ALT_NON_REF = 5;
const BYTE REF_ALT_N = 6;


//                                     A  T  G  C  N
const BYTE TRANSITION_MAP_BASE_A[5] = {0, 0, 1, 0, 0}; // A->G
const BYTE TRANSITION_MAP_BASE_T[5] = {0, 0, 0, 1, 0}; // T->C
const BYTE TRANSITION_MAP_BASE_G[5] = {1, 0, 0, 0, 0}; // G->A
const BYTE TRANSITION_MAP_BASE_C[5] = {0, 1, 0, 0, 0}; // C->T
const BYTE TRANSITION_MAP_BASE_N[5] = {0, 0, 0, 0, 0}; // None
const BYTE* const TRANSITION_MAP[9] =
		{TRANSITION_MAP_BASE_A,
		 TRANSITION_MAP_BASE_T,
		 TRANSITION_MAP_BASE_G,
		 TRANSITION_MAP_BASE_C,
		 TRANSITION_MAP_BASE_N,
		 TRANSITION_MAP_BASE_N,
		 TRANSITION_MAP_BASE_N,
		 TRANSITION_MAP_BASE_N,
		 TRANSITION_MAP_BASE_N};

//                                       A  T  G  C  N
const BYTE TRANSVERSION_MAP_BASE_A[5] = {0, 1, 0, 1, 0};
const BYTE TRANSVERSION_MAP_BASE_T[5] = {1, 0, 1, 0, 0};
const BYTE TRANSVERSION_MAP_BASE_G[5] = {0, 1, 0, 1, 0};
const BYTE TRANSVERSION_MAP_BASE_C[5] = {1, 0, 1, 0, 0};
const BYTE TRANSVERSION_MAP_BASE_N[5] = {0, 0, 0, 0, 0};
const BYTE* const TRANSVERSION_MAP[9] =
		{TRANSVERSION_MAP_BASE_A,
		 TRANSVERSION_MAP_BASE_T,
		 TRANSVERSION_MAP_BASE_G,
		 TRANSVERSION_MAP_BASE_C,
		 TRANSVERSION_MAP_BASE_N,
		 TRANSVERSION_MAP_BASE_N,
		 TRANSVERSION_MAP_BASE_N,
		 TRANSVERSION_MAP_BASE_N,
		 TRANSVERSION_MAP_BASE_N};

/*------   EOF markers   ------*/
// EOF markers
// 64-bit XXHASH of "TACHYON-BLOCK-EOF"
const U64 TACHYON_BLOCK_EOF = 7964708207515128046;
// 32b SHA-256 digest
// echo -n "TACHYON-EOF" | openssl dgst -sha256
const std::string TACHYON_FILE_EOF = "fa20427e118a1f0c7e1195e4f85e15e74368b112e70dd89fd02772e1d90cb5dd";
const BYTE TACHYON_FILE_EOF_LENGTH = 32;

}
}

#endif /* MAGICCONSTANTS_H_ */
