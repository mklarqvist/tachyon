#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "containers/variant_block.h"
#include "containers/stride_container.h"

namespace tachyon{
namespace algorithm{

class CompressionManager{
public:
	typedef CompressionManager        self_type;
	typedef UncompressedCodec         no_codec_type;
	typedef ZSTDCodec                 zstd_codec_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer container_type;
	typedef containers::VariantBlockFooter footer_type;

public:
	CompressionManager() = default;
	~CompressionManager() = default;

	bool Compress(variant_block_type& block, const uint8_t general_level = 6);
	bool Decompress(variant_block_type& block);
	bool Decompress(container_type& container, yon_gt_ppa& gt_ppa);

	template <class T>
	bool StrideTranspose(container_type& container){
		if(container.header.data_header.IsUniform() == false &&
		   (container.header.data_header.GetPrimitiveType() == YON_TYPE_32B || container.header.data_header.GetPrimitiveType() == YON_TYPE_16B || container.header.data_header.GetPrimitiveType() == YON_TYPE_8B))
		{
			io::BasicBuffer temp(container.data_uncompressed.size() + 65536);
			containers::StrideContainer<> s(container);
			const T* data = reinterpret_cast<const T*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(T);
			assert(container.data_uncompressed.size() % sizeof(T) == 0);
			assert(n_entries == s.size());

			for(int i = 0; i < s.size(); ++i){

			}
		}
		return true;
	}


	/**<
	 * Decompress a data container.
	 * @param container Src and dst target container.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool Decompress(container_type& container);

public:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
