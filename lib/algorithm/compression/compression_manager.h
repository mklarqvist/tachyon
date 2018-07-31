#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "containers/variant_block.h"

namespace tachyon{
namespace algorithm{

class CompressionManager{
private:
	typedef CompressionManager        self_type;
	typedef UncompressedCodec         no_codec_type;
	typedef ZSTDCodec                 zstd_codec_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer container_type;

public:
	CompressionManager() = default;
	~CompressionManager() = default;

	bool Compress(variant_block_type& block, const BYTE general_level = 6, const BYTE float_level = 3);
	bool Decompress(variant_block_type& block);
	bool Decompress(algorithm::PermutationManager& permutation_manager);

	/**<
	 * Decompress an abstract data container
	 * @param container Target container
	 * @return          Returns TRUE upon success or FALSE otherwise
	 */
	bool Decompress(container_type& container);

public:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
