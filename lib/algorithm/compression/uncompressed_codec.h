#ifndef ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_
#define ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_

#include "compression_container.h"

namespace tachyon{
namespace algorithm{

class UncompressedCodec : public CompressionContainer{
public:
	typedef UncompressedCodec self_type;
	typedef yon1_dc_t    container_type;
	typedef yon_buffer_t buffer_type;
	typedef yon_gt_ppa   permutation_type;

public:
	UncompressedCodec() = default;
	~UncompressedCodec() = default;

	bool Compress(container_type& container, permutation_type& manager) { return true; }
	bool Compress(container_type& container);
	inline bool CompressStrides(container_type& container) { return true; }
	bool Decompress(container_type& container);
	bool DecompressStrides(container_type& container);

protected:
	buffer_type buffer;
};

}
}



#endif /* ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_ */
