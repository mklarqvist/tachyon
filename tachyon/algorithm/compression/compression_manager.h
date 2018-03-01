#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "../../containers/variantblock.h"
#include "compression_container.h"

namespace tachyon{
namespace algorithm{

class CompressionManager{
private:
	typedef CompressionManager        self_type;
	typedef UncompressedCodec         no_codec_type;
	typedef ZSTDCodec                 zstd_codec_type;
	typedef containers::VariantBlock  block_type;
	typedef containers::DataContainer container_type;

public:
	CompressionManager(){}
	~CompressionManager(){}

	bool compress(block_type& block){
		zstd_codec.setCompressionLevel(20);
		zstd_codec.setCompressionLevelData(3);
		if(block.index_entry.controller.hasGTPermuted) zstd_codec.encode(block.ppa_manager);

		zstd_codec.setCompressionLevel(20);
		if(block.meta_hot_container.n_entries)        zstd_codec.encode(block.meta_hot_container);
		if(block.gt_rle_container.n_entries)          zstd_codec.encode(block.gt_rle_container);
		if(block.gt_simple_container.n_entries)       zstd_codec.encode(block.gt_simple_container);
		if(block.gt_support_data_container.n_entries) zstd_codec.encode(block.gt_support_data_container);
		if(block.meta_cold_container.n_entries)       zstd_codec.encode(block.meta_cold_container);
		if(block.meta_info_map_ids.n_entries)         zstd_codec.encode(block.meta_info_map_ids);
		if(block.meta_filter_map_ids.n_entries)       zstd_codec.encode(block.meta_filter_map_ids);
		if(block.meta_format_map_ids.n_entries)       zstd_codec.encode(block.meta_format_map_ids);

		for(U32 i = 0; i < block.index_entry.n_info_streams; ++i){
			if(block.info_containers[i].header.controller.type == YON_TYPE_FLOAT ||
			   block.info_containers[i].header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(20);
			}
			else zstd_codec.setCompressionLevel(20);
			zstd_codec.encode(block.info_containers[i]);
		}

		for(U32 i = 0; i < block.index_entry.n_format_streams; ++i){
			if(block.format_containers[i].header.controller.type == YON_TYPE_FLOAT ||
			   block.format_containers[i].header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(20);
			}
			else zstd_codec.setCompressionLevel(20);
			zstd_codec.encode(block.format_containers[i]);
		}

		return true;
	}

	bool decompress(block_type& block){
		if(block.ppa_manager.PPA.size()){
			if(!this->decompress(block.ppa_manager)){
				std::cerr << "failed to decompress ppa!" << std::endl;
			}
		}

		if(block.meta_hot_container.getSizeCompressed()){
			//std::cerr << "decompress hot: " << block.meta_hot_container.buffer_data.size() << std::endl;
			if(!this->decompress(block.meta_hot_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_cold_container.getSizeCompressed()){
			//std::cerr << "decompress cold" << std::endl;
			if(!this->decompress(block.meta_cold_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_rle_container.getSizeCompressed()){
			//std::cerr << "decompress gt rle" << std::endl;
			if(!this->decompress(block.gt_rle_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_simple_container.getSizeCompressed()){
			//std::cerr << "decompress gt simple" << std::endl;
			if(!this->decompress(block.gt_simple_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.gt_support_data_container.getSizeCompressed()){
			//std::cerr << "decompress gt support" << std::endl;
			if(!this->decompress(block.gt_support_data_container)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_info_map_ids.getSizeCompressed()){
			//std::cerr << "decompress info map ids" << std::endl;
			if(!this->decompress(block.meta_info_map_ids)){
				std::cerr << "failed to meta_info_map_ids!" << std::endl;
			}
		}

		if(block.meta_filter_map_ids.getSizeCompressed()){
			//std::cerr << "decompress filter map ids" << std::endl;
			if(!this->decompress(block.meta_filter_map_ids)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		if(block.meta_format_map_ids.getSizeCompressed()){
			//std::cerr << "decompress format map ids" << std::endl;
			if(!this->decompress(block.meta_format_map_ids)){
				std::cerr << "failed to decompress!" << std::endl;
			}
		}

		for(U32 i = 0; i < block.index_entry.n_info_streams; ++i){
			//std::cerr << "decompress: " << i << "/" << block.index_entry.n_info_streams << std::endl;
			if(block.info_containers[i].getSizeCompressed()){
				if(!this->decompress(block.info_containers[i])){
					std::cerr << "failed to decompress!" << std::endl;
				}
			}
		}

		for(U32 i = 0; i < block.index_entry.n_format_streams; ++i){
			//std::cerr << "decompress: " << i << "/" << block.index_entry.n_format_streams << std::endl;
			if(block.format_containers[i].getSizeCompressed()){
				if(!this->decompress(block.format_containers[i])){
					std::cerr << "failed to decompress!" << std::endl;
				}
			}
		}

		//std::cerr << "returning from decompressing" << std::endl;
		return true;
	}

	bool decompress(algorithm::PermutationManager& permutation_manager){
		return(this->zstd_codec.decode(permutation_manager));
	}

	bool decompress(container_type& container){
		if(container.header.controller.encoder == YON_ENCODE_ZSTD){
			if(!this->zstd_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else if(container.header.controller.encoder == YON_ENCODE_NONE){
			if(!this->no_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else {
			std::cerr << "ILLEGAL CODEC" << std::endl;
			return false;
		}

		if(container.header.controller.mixedStride){
			if(container.header_stride.controller.encoder == YON_ENCODE_ZSTD){
				this->zstd_codec.decodeStrides(container);
			} else if (container.header_stride.controller.encoder == YON_ENCODE_NONE){
				this->no_codec.decodeStrides(container);
			}
		}
		return true;
	}

private:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
