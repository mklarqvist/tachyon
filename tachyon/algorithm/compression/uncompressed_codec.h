#ifndef ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_
#define ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_

#include "compression_container.h"

namespace tachyon{
namespace algorithm{

class UncompressedCodec : public CompressionContainer{
private:
	typedef UncompressedCodec             self_type;

protected:
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	UncompressedCodec() = default;
	~UncompressedCodec() = default;

	inline const bool compress(permutation_type& manager){ return true; }
	const bool compress(container_type& container);
	inline const bool compressStrides(container_type& container){ return true; }
	const bool decompress(container_type& container);
	const bool decompressStrides(container_type& container);

protected:
	buffer_type buffer;
};

}
}



#endif /* ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_ */
