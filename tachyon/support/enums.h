#ifndef SUPPORT_ENUMS_H_
#define SUPPORT_ENUMS_H_

namespace tachyon{

enum TACHYON_CORE_TYPE{
	YON_TYPE_UNKNOWN,
	YON_TYPE_8B,
	YON_TYPE_16B,
	YON_TYPE_32B,
	YON_TYPE_64B,
	YON_TYPE_FLOAT,
	YON_TYPE_DOUBLE,
	YON_TYPE_BOOLEAN,
	YON_TYPE_CHAR,
	YON_TYPE_STRUCT
};

enum TACHYON_CORE_COMPRESSION{
	YON_ENCODE_NONE,
	YON_ENCODE_ZSTD,
	YON_ENCODE_ZPAQ
};


enum TACHYON_ENCRYPTION{YON_ENCRYPTION_NONE,
	                    YON_ENCRYPTION_AES_128,
					    YON_ENCRYPTION_AES_256_GCM,
					    YON_ENCRYPTION_RSA_4096};

}

/**<
 * For genotype encodings: describes what encoding
 * algorithm was used to store genotypes
 */
enum TACHYON_GT_ENCODING{
	YON_GT_RLE_DIPLOID_BIALLELIC = 0,//!< YON_GT_RLE_DIPLOID_BIALLELIC
	YON_GT_RLE_DIPLOID_NALLELIC  = 1,//!< YON_GT_RLE_DIPLOID_NALLELIC
	YON_GT_BCF_DIPLOID           = 2,//!< YON_GT_BCF_DIPLOID
	YON_GT_BCF_STYLE             = 3 //!< YON_GT_BCF_STYLE
};

/**<
 * For genotype encodings: this enum describes what
 * primitive type the data is stored as
 */
enum TACHYON_GT_PRIMITIVE_TYPE{
	YON_GT_BYTE = 0,//!< YON_GT_BYTE
	YON_GT_U16  = 1,//!< YON_GT_U16
	YON_GT_U32  = 2,//!< YON_GT_U32
	YON_GT_U64  = 3 //!< YON_GT_U64
};



#endif /* SUPPORT_ENUMS_H_ */
