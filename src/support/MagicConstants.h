#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>
#include "../support/TypeDefinitions.h"

extern int SILENT;

namespace Tachyon{
namespace Constants{

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

// Versioning
const float PROGRAM_VERSION = 0.4; // major

const double ALLOWED_ROUNDING_ERROR = 0.001;

const std::string PROGRAM_NAME = "tachyon";
const std::string OUTPUT_SUFFIX = "tyn";

// Headers
const char* const WRITE_HEADER_MAGIC = "TOTEMPOLE\1";
const U16 WRITE_HEADER_MAGIC_LENGTH = 10;

const BYTE TOMAHAWK_ALLELE_PACK_WIDTH = 2; // bit / allele
const BYTE TOMAHAWK_SNP_PACK_WIDTH = TOMAHAWK_ALLELE_PACK_WIDTH * 2; // bits / genotype
const BYTE TOMAHAWK_SHIFT_SIZE = Constants::TOMAHAWK_SNP_PACK_WIDTH + 1;

// Encoding for alleles
const char TOMAHAWK_ALLELE_LOOKUP[4] = {2, 3, 0, 1};
const char TOMAHAWK_ALLELE_LOOKUP_REVERSE[4] = {'0', '1', '.', '?'};

// 0 -> 0
// 1 -> 1
// 2 -> 4
// 3 -> 5
const BYTE TOMAHAWK_ALLELE_REDUCED_MAP[4] = {0, 1, 4, 5};
// 0 -> 0
// 1 -> 1
// 4 -> 4
// 5 -> 5
const BYTE TOMAHAWK_ALLELE_SELF_MAP[6] = {0, 1, 0, 0, 4, 5};

// Encoding for bases
const char* const REF_ALT_LOOKUP = "ATGCN";
const BYTE REF_ALT_A = 0;
const BYTE REF_ALT_T = 1;
const BYTE REF_ALT_G = 2;
const BYTE REF_ALT_C = 3;
const BYTE REF_ALT_N = 4;

//EOF
const BYTE eof_length = 6;
const U64 eof[6] = {2336361506924422487, 7959953386435011938, 8243124871055238688, 2334386829831791136, 8583987794834190964, 28464622577219173};
// EOF poem: "We will be known forever by the tracks we leave"

}
}

#endif /* MAGICCONSTANTS_H_ */
