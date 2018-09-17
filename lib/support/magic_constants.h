#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>
#include <regex>

extern int SILENT;

namespace tachyon{

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

/*------   Version   ------*/
const int32_t TACHYON_VERSION_MAJOR = 0;
const int32_t TACHYON_VERSION_MINOR = 5;
const int32_t TACHYON_VERSION_PATCH = 0;
const int32_t TACHYON_VERSION_NUMBER  = (TACHYON_VERSION_MAJOR *100*100 + TACHYON_VERSION_MINOR *100 + TACHYON_VERSION_PATCH);
const std::string TACHYON_LIB_VERSION = std::to_string(TACHYON_VERSION_MAJOR) + '.' + std::to_string(TACHYON_VERSION_MINOR) + '.' + std::to_string(TACHYON_VERSION_PATCH);

/*------   Basics   ------*/
const std::string TACHYON_PROGRAM_NAME  = "tachyon";
const std::string TACHYON_OUTPUT_SUFFIX = "yon";
const std::string TACHYON_MAGIC_HEADER  = "TACHYON\1";
const uint32_t    TACHYON_MAGIC_HEADER_LENGTH = 8;

/*------   Regular expression patterns  ------*/
const std::regex YON_REGEX_CANONICAL_BASES = std::regex("^[ATGC]+$");
const std::regex YON_REGEX_CONTIG_ONLY     = std::regex("^[A-Za-z0-9\\-_]+$");
const std::regex YON_REGEX_CONTIG_POSITION = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?$");
const std::regex YON_REGEX_CONTIG_RANGE    = std::regex("^[A-Za-z0-9\\-_]+\\:[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?\\-[0-9]+([\\.]{1}[0-9]+){0,1}([eE]{1}[0-9]{1})?$");
const std::regex YON_REGEX_PACKED_ALLELES  = std::regex("^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$");

/*------   EOF markers   ------*/
// 64-bit XXHASH of "TACHYON-BLOCK-EOF"
const uint64_t TACHYON_BLOCK_EOF = 7964708207515128046;
// 32b SHA-256 digest
// echo -n "TACHYON-EOF" | openssl dgst -sha256
const std::string TACHYON_FILE_EOF = "fa20427e118a1f0c7e1195e4f85e15e74368b112e70dd89fd02772e1d90cb5dd";
const uint32_t TACHYON_FILE_EOF_LENGTH = 32;

}

#endif /* MAGICCONSTANTS_H_ */
