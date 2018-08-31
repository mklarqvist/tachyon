#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "containers/variant_block.h"

// Hide logic. Doesn't compile otherwise.
namespace FastPForLib{
class IntegerCODEC;
}

namespace tachyon{
namespace algorithm{

struct yon_comp_memory_pair {
	yon_comp_memory_pair() : codec(0), times_seen(10){}

	uint8_t codec;
	uint32_t times_seen;
};

class CompressionManager{
private:
	typedef CompressionManager        self_type;
	typedef UncompressedCodec         no_codec_type;
	typedef ZSTDCodec                 zstd_codec_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer container_type;

public:
	CompressionManager();
	~CompressionManager();

	bool Compress(variant_block_type& block, const uint8_t general_level = 6, const uint8_t float_level = 3);
	bool Decompress(variant_block_type& block);
	bool Decompress(container_type& container, yon_gt_ppa& gt_ppa);

	int32_t CompressCodecWrapper(FastPForLib::IntegerCODEC* codec, container_type& container, const bool move = true);
	int32_t CompressCodec1(container_type& container, const bool move = true);
	int32_t CompressCodec2(container_type& container, const bool move = true);
	int32_t CompressCodec3(container_type& container, const bool move = true);
	int32_t CompressCodec4(container_type& container, const bool move = true);
	int32_t CompressDelta(container_type& container, const bool move = true);
	bool CompressEvaluate(container_type& container, yon_comp_memory_pair* memory, const uint32_t global_id);

	bool CompressEvaluate(container_type& container, yon_gt_ppa* ppa);


	/**<
	 * Decompress an abstract data container
	 * @param container Target container
	 * @return          Returns TRUE upon success or FALSE otherwise
	 */
	bool Decompress(container_type& container);

public:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
	yon_comp_memory_pair memory_basic[100];
	yon_comp_memory_pair memory_info[100];
	yon_comp_memory_pair memory_format[100];

	io::BasicBuffer backup_buffer;
	io::BasicBuffer support_buffer;
	FastPForLib::IntegerCODEC* codec1;
	FastPForLib::IntegerCODEC* codec2;
	FastPForLib::IntegerCODEC* codec3;
	FastPForLib::IntegerCODEC* codec4;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
