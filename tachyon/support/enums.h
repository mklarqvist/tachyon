#ifndef SUPPORT_ENUMS_H_
#define SUPPORT_ENUMS_H_

namespace tachyon{

enum TACHYON_CORE_TYPE{
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
	YON_ENCODE_ZSTD
};

}



#endif /* SUPPORT_ENUMS_H_ */
