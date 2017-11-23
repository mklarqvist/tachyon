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
const float PROGRAM_VERSION = 0.1; // major
const std::string PROGRAM_NAME = "tachyon";
const std::string OUTPUT_SUFFIX = "tyn";

const BYTE ALLELE_PACK_WIDTH = 2; // bit / allele
const BYTE SNP_PACK_WIDTH = ALLELE_PACK_WIDTH * 2; // bits / genotype
const BYTE SHIFT_SIZE = Constants::SNP_PACK_WIDTH + 1;

// Encoding for alleles
const char ALLELE_LOOKUP[4] = {2, 3, 0, 1};
const char ALLELE_LOOKUP_REVERSE[4] = {'0', '1', '.', '?'};

// 0 -> 0
// 1 -> 1
// 2 -> 4
// 3 -> 5
const BYTE ALLELE_REDUCED_MAP[4] = {0, 1, 4, 5};
// 0 -> 0
// 1 -> 1
// 4 -> 4
// 5 -> 5
const BYTE ALLELE_SELF_MAP[6] = {0, 1, 0, 0, 4, 5};

// Encoding for bases
const char* const REF_ALT_LOOKUP = "ATGCN";
const BYTE REF_ALT_A = 0;
const BYTE REF_ALT_T = 1;
const BYTE REF_ALT_G = 2;
const BYTE REF_ALT_C = 3;
const BYTE REF_ALT_N = 4;

// BLOCK EOF
const U64 BLOCK_EOF = 7964708207515128046;
const U64 FILE_EOF = 1404072766443001870;

}
}

#endif /* MAGICCONSTANTS_H_ */
