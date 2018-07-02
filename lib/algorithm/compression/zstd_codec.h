#ifndef ALGORITHM_COMPRESSION_ZSTD_CODEC_H_
#define ALGORITHM_COMPRESSION_ZSTD_CODEC_H_

#include "zstd.h"
#include "zstd_errors.h"
#include "compression_container.h"

namespace tachyon{
namespace algorithm{

class ZSTDCodec : public CompressionContainer{
private:
	typedef        ZSTDCodec   self_type;
	typedef struct ZSTD_CCtx_s ZSTD_CCtx;
	typedef struct ZSTD_DCtx_s ZSTD_DCtx;

public:
	ZSTDCodec();
	~ZSTDCodec();

	inline void setCompressionLevel(const S32& c){ this->compression_level_data = c; this->compression_level_strides = c; }
	inline void setCompressionLevelData(const S32& c){ this->compression_level_data = c; }
	inline void setCompressionLevelStrides(const S32& c){ this->compression_level_strides = c; }

	const bool compress(container_type& container);
	const bool compressStrides(container_type& container);
	const bool compress(permutation_type& manager);
	const bool decompress(container_type& container);
	const bool decompressStrides(container_type& container);
	const bool decompress(permutation_type& manager);

private:
	S32 compression_level_data;
	S32 compression_level_strides;
	ZSTD_CCtx* compression_context_;
	ZSTD_DCtx* decompression_context_;
};

}
}



#endif /* ALGORITHM_COMPRESSION_ZSTD_CODEC_H_ */
