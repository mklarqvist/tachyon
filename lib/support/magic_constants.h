#ifndef MAGICCONSTANTS_H_
#define MAGICCONSTANTS_H_

#include <string>
#include <regex>

extern int SILENT;

/**<
 * These definitions correspond to the array offsets for the
 * invariant containers in the yon1_vb_t/yon_vb_ftr.
 * These correspond to the core components for site information
 * and internals such as CONTROLLER, ID_*, and GT_SUPPORT and
 * GT_PLOIDY.
 */
#define YON_BLK_N_STATIC   25// Total number of invariant headers
#define YON_BLK_PPA         0 // Sample permutation array
#define YON_BLK_CONTIG      1
#define YON_BLK_POSITION    2
#define YON_BLK_REFALT      3
#define YON_BLK_CONTROLLER  4 // Set memberships
#define YON_BLK_QUALITY     5
#define YON_BLK_NAMES       6
#define YON_BLK_ALLELES     7
#define YON_BLK_ID_INFO     8
#define YON_BLK_ID_FORMAT   9
#define YON_BLK_ID_FILTER  10
#define YON_BLK_GT_INT8    11 // Run-length encoded genotypes
#define YON_BLK_GT_INT16   12
#define YON_BLK_GT_INT32   13
#define YON_BLK_GT_INT64   14
#define YON_BLK_GT_S_INT8  15 // Standard encoded genotypes
#define YON_BLK_GT_S_INT16 16
#define YON_BLK_GT_S_INT32 17
#define YON_BLK_GT_S_INT64 18
#define YON_BLK_GT_N_INT8  19 // Standard encoded genotypes
#define YON_BLK_GT_N_INT16 20
#define YON_BLK_GT_N_INT32 21
#define YON_BLK_GT_N_INT64 22
#define YON_BLK_GT_SUPPORT 23 // Genotype support
#define YON_BLK_GT_PLOIDY  24 // Genotype ploidy

#define YON_BLK_BV_PPA         (1 << YON_BLK_PPA)
#define YON_BLK_BV_CONTIG      (1 << YON_BLK_CONTIG)
#define YON_BLK_BV_POSITION    (1 << YON_BLK_POSITION)
#define YON_BLK_BV_REFALT      (1 << YON_BLK_REFALT)
#define YON_BLK_BV_CONTROLLER  (1 << YON_BLK_CONTROLLER)
#define YON_BLK_BV_QUALITY     (1 << YON_BLK_QUALITY)
#define YON_BLK_BV_NAMES       (1 << YON_BLK_NAMES)
#define YON_BLK_BV_ALLELES     (1 << YON_BLK_ALLELES)
#define YON_BLK_BV_ID_INFO     (1 << YON_BLK_ID_INFO)
#define YON_BLK_BV_ID_FORMAT   (1 << YON_BLK_ID_FORMAT)
#define YON_BLK_BV_ID_FILTER   (1 << YON_BLK_ID_FILTER)
#define YON_BLK_BV_GT_INT8     (1 << YON_BLK_GT_INT8)
#define YON_BLK_BV_GT_INT16    (1 << YON_BLK_GT_INT16)
#define YON_BLK_BV_GT_INT32    (1 << YON_BLK_GT_INT32)
#define YON_BLK_BV_GT_INT64    (1 << YON_BLK_GT_INT64)
#define YON_BLK_BV_GT_S_INT8   (1 << YON_BLK_GT_S_INT8)
#define YON_BLK_BV_GT_S_INT16  (1 << YON_BLK_GT_S_INT16)
#define YON_BLK_BV_GT_S_INT32  (1 << YON_BLK_GT_S_INT32)
#define YON_BLK_BV_GT_S_INT64  (1 << YON_BLK_GT_S_INT64)
#define YON_BLK_BV_GT_N_INT8   (1 << YON_BLK_GT_N_INT8)
#define YON_BLK_BV_GT_N_INT16  (1 << YON_BLK_GT_N_INT16)
#define YON_BLK_BV_GT_N_INT32  (1 << YON_BLK_GT_N_INT32)
#define YON_BLK_BV_GT_N_INT64  (1 << YON_BLK_GT_N_INT64)
#define YON_BLK_BV_GT_SUPPORT  (1 << YON_BLK_GT_SUPPORT)
#define YON_BLK_BV_GT_PLOIDY   (1 << YON_BLK_GT_PLOIDY)
// Special values
#define YON_BLK_BV_INFO        (1 << (YON_BLK_N_STATIC))
#define YON_BLK_BV_FORMAT      (1 << (YON_BLK_N_STATIC + 1))
#define YON_BLK_BV_GT          ((YON_BLK_BV_GT_INT8)|(YON_BLK_BV_GT_INT16)|(YON_BLK_BV_GT_INT32)|(YON_BLK_BV_GT_INT64)|(YON_BLK_BV_GT_S_INT8)|(YON_BLK_BV_GT_S_INT16)|(YON_BLK_BV_GT_S_INT32)|(YON_BLK_BV_GT_S_INT64)|(YON_BLK_BV_GT_N_INT8)|(YON_BLK_BV_GT_N_INT16)|(YON_BLK_BV_GT_N_INT32)|(YON_BLK_BV_GT_N_INT64)|(YON_BLK_BV_GT_SUPPORT)|(YON_BLK_BV_GT_PLOIDY))

namespace tachyon {

const std::vector<std::string> YON_BLK_PRINT_NAMES = {
		"PPA","CONTIG","POSITION","REFALT","CONTROLLER","QUALITY","NAMES",
		"ALLELES","ID_INFO","ID_FORMAT","ID_FILTER",
		"GT_INT8","GT_INT16","GT_INT32","GT_INT64",
		"GT_S_INT8","GT_S_INT16","GT_S_INT32","GT_S_INT64",
		"GT_N_INT8","GT_N_INT16","GT_N_INT32","GT_N_INT64",
		"GT_SUPPORT","GT_PLOIDY"};

extern std::string LITERAL_COMMAND_LINE;
extern std::string INTERPRETED_COMMAND;

/*------   Version   ------*/
const int32_t TACHYON_VERSION_MAJOR = 0;
const int32_t TACHYON_VERSION_MINOR = 6;
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

/*------ Core enums --------*/
typedef enum {
	YON_TYPE_UNKNOWN, YON_TYPE_8B, YON_TYPE_16B, YON_TYPE_32B, YON_TYPE_64B,
	YON_TYPE_FLOAT, YON_TYPE_DOUBLE, YON_TYPE_BOOLEAN, YON_TYPE_CHAR,
	YON_TYPE_STRUCT
} TACHYON_CORE_TYPE;

typedef enum { YON_ENCODE_NONE, YON_ENCODE_ZSTD, YON_ENCODE_ZPAQ } TACHYON_CORE_COMPRESSION;

typedef enum {
	YON_ENCRYPTION_NONE, YON_ENCRYPTION_AES_128, YON_ENCRYPTION_AES_256_GCM,
	YON_ENCRYPTION_RSA_4096
} TACHYON_ENCRYPTION;

typedef enum {
	YON_GT_RLE_DIPLOID_BIALLELIC, YON_GT_RLE_DIPLOID_NALLELIC, YON_GT_BCF_DIPLOID,
	YON_GT_BCF_STYLE, YON_GT_RLE_NPLOID
} TACHYON_GT_ENCODING;

typedef enum { YON_GT_BYTE, YON_GT_U16, YON_GT_U32, YON_GT_U64 } TACHYON_GT_PRIMITIVE_TYPE;

typedef enum {
	YON_VCF_HEADER_INTEGER, YON_VCF_HEADER_FLOAT, YON_VCF_HEADER_FLAG,
	YON_VCF_HEADER_CHARACTER, YON_VCF_HEADER_STRING
} TACHYON_VARIANT_HEADER_FIELD_TYPE;

typedef enum {
	YON_VARIANT_CLASS_SNP, YON_VARIANT_CLASS_MNP, YON_VARIANT_CLASS_INDEL,
	YON_VARIANT_CLASS_SV, YON_VARIANT_CLASS_CLUMPED, YON_VARIANT_CLASS_UNKNOWN
} TACHYON_VARIANT_CLASSIFICATION_TYPE;

const std::string TACHYON_VARIANT_CLASSIFICATION_STRING[6] = {"SNP","MNP","INDEL","SV","CLUMPED","UNKNOWN"};

typedef enum {
	YON_CMP_LESS, YON_CMP_LESS_EQUAL, YON_CMP_GREATER, YON_CMP_GREATER_EQUAL,
	YON_CMP_EQUAL, YON_CMP_NOT_EQUAL, YON_CMP_REGEX
} TACHYON_COMPARATOR_TYPE;

typedef enum {
	YON_FILTER_NUMBER_ALT_ALLELES,
	YON_FILTER_MIXED_PHASING,
	YON_FILTER_MIXED_PLOIDY,
	YON_FILTER_MISSING_GT,
	YON_FILTER_ALLELE_FREQUENCY,
	YON_FILTER_ALLELE_COUNT,
	YON_FILTER_UNIFORM_PHASE,
	YON_FILTER_KNOWN_NOVEL,
	YON_FILTER_REFERENCE_ALLELE,
	YON_FILTER_ALT_ALLELE,
	YON_FILTER_NAME,
	YON_FILTER_UNSEEN_ALT,
	YON_FILTER_QUALITY
} TACHYON_FILTER_FUNCTION;

}

#endif /* MAGICCONSTANTS_H_ */
