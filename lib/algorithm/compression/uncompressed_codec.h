#ifndef ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_
#define ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_

#include "compression_container.h"

namespace tachyon{
namespace algorithm{

class UncompressedCodec : public CompressionContainer{
private:
	typedef UncompressedCodec             self_type;
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	UncompressedCodec() = default;
	~UncompressedCodec() = default;

	inline bool compress(permutation_type& manager){ return true; }
	bool compress(container_type& container);
	inline bool compressStrides(container_type& container){ return true; }
	bool decompress(container_type& container);
	bool decompressStrides(container_type& container);

protected:
	buffer_type buffer;
};

}
}



#endif /* ALGORITHM_COMPRESSION_UNCOMPRESSED_CODEC_H_ */
