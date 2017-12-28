#ifndef ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_

#include "CompressionContainer.h"
#include "../../core/StreamContainer.h"

namespace Tachyon{
namespace Compression{

class CompressionManager{
private:
	typedef CompressionManager self_type;
	typedef UncompressedCodec no_codec_type;
	typedef ZSTDCodec zstd_codec_type;

	typedef Core::BlockEntry block_type;
	typedef Core::StreamContainer container_type;

public:
	CompressionManager(){}
	~CompressionManager(){}

	bool decompress(block_type& block){
		if(block.meta_hot_container.sizeCompressed()){
			if(!this->decompress(block.meta_hot_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_cold_container.sizeCompressed()){
			if(!this->decompress(block.meta_cold_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_rle_container.sizeCompressed()){
			if(!this->decompress(block.gt_rle_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_simple_container.sizeCompressed()){
			if(!this->decompress(block.gt_simple_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_support_data_container.sizeCompressed()){
			if(!this->decompress(block.gt_support_data_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_info_map_ids.sizeCompressed()){
			if(!this->decompress(block.meta_info_map_ids)){
				std::cerr << "failed to meta_info_map_ids!" << std::endl;
			}
		}

		if(block.meta_filter_map_ids.sizeCompressed()){
			if(!this->decompress(block.meta_filter_map_ids)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_format_map_ids.sizeCompressed()){
			if(!this->decompress(block.meta_format_map_ids)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		for(U32 i = 0; i < block.index_entry.n_info_streams; ++i){
			if(block.info_containers[i].sizeCompressed()){
				if(!this->decompress(block.info_containers[i])){
					std::cerr << "failed to decompress!" << std::endl;
				}
			}
		}

		for(U32 i = 0; i < block.index_entry.n_format_streams; ++i){
			if(block.format_containers[i].sizeCompressed()){
				if(!this->decompress(block.format_containers[i])){
					std::cerr << "failed to decompress!" << std::endl;
				}
			}
		}

		return true;
	}

	bool decompress(container_type& container){
		if(container.header.controller.encoder == Core::ENCODE_ZSTD){
			if(!this->zstd_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else if(container.header.controller.encoder == Core::ENCODE_NONE){
			if(!this->no_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else {
			std::cerr << "ILLEGAL CODEC" << std::endl;
			return false;
		}

		if(container.header.controller.mixedStride){
			if(container.header_stride.controller.encoder == Core::ENCODE_ZSTD){
				this->zstd_codec.decodeStrides(container);
			} else if (container.header_stride.controller.encoder == Core::ENCODE_NONE){
				this->no_codec.decodeStrides(container);
			}
		}
		return true;
	}

private:
	no_codec_type no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_ */
