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

	inline void SetCompressionLevel(const S32& c){ this->compression_level_data = c; this->compression_level_strides = c; }
	inline void SetCompressionLevelData(const S32& c){ this->compression_level_data = c; }
	inline void SetCompressionLevelStrides(const S32& c){ this->compression_level_strides = c; }

	bool Compress(container_type& container);
	bool CompressStrides(container_type& container);
	bool Compress(permutation_type& manager);
	bool Decompress(container_type& container);
	bool DecompressStrides(container_type& container);
	bool Decompress(permutation_type& manager);

	bool Compress(const io::BasicBuffer& src, io::BasicBuffer& dst, const int compression_level);
	bool Decompress(const io::BasicBuffer& src, io::BasicBuffer& dst);

private:
	S32 compression_level_data;
	S32 compression_level_strides;
	ZSTD_CCtx* compression_context_;
	ZSTD_DCtx* decompression_context_;
};

}
}



#endif /* ALGORITHM_COMPRESSION_ZSTD_CODEC_H_ */
