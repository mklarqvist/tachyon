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
	typedef ZPAQContainer             zpaq_codec_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer container_type;

public:
	CompressionManager(){}
	~CompressionManager(){}

	bool compress(variant_block_type& block){
		zstd_codec.setCompressionLevel(20);
		zstd_codec.setCompressionLevelData(3);
		if(block.header.controller.hasGTPermuted) zstd_codec.compress(block.ppa_manager);

		zstd_codec.setCompressionLevel(20);
		if(block.meta_hot_container.n_entries)        zstd_codec.compress(block.meta_hot_container);
		if(block.gt_rle_container.n_entries)          zstd_codec.compress(block.gt_rle_container);
		if(block.gt_simple_container.n_entries)       zstd_codec.compress(block.gt_simple_container);
		if(block.gt_support_data_container.n_entries) zstd_codec.compress(block.gt_support_data_container);
		if(block.meta_cold_container.n_entries)       zstd_codec.compress(block.meta_cold_container);
		if(block.meta_info_map_ids.n_entries)         zstd_codec.compress(block.meta_info_map_ids);
		if(block.meta_filter_map_ids.n_entries)       zstd_codec.compress(block.meta_filter_map_ids);
		if(block.meta_format_map_ids.n_entries)       zstd_codec.compress(block.meta_format_map_ids);

		for(U32 i = 0; i < block.header.n_info_streams; ++i){
			if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(20);
				//zpaq_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
			}
			else {
				zstd_codec.setCompressionLevel(20);
				//zpaq_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
			}
			zstd_codec.compress(block.info_containers[i]);
			//zpaq_codec.compress(block.info_containers[i], "4");
		}

		for(U32 i = 0; i < block.header.n_format_streams; ++i){
			if(block.format_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.format_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(20);
				//zpaq_codec.compress(block.format_containers[i]);
				//zstd_codec.compress(block.format_containers[i]);
			}
			else {
				zstd_codec.setCompressionLevel(20);
				//zstd_codec.compress(block.format_containers[i]);
				//zpaq_codec.compress(block.format_containers[i]);
			}
			zstd_codec.compress(block.format_containers[i]);
		}

		return true;
	}

	bool decompress(variant_block_type& block){
		if(block.ppa_manager.PPA.size()){
			if(!this->decompress(block.ppa_manager)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT permutation information!" << std::endl;
			}
		}

		if(block.meta_hot_container.getSizeCompressed()){
			if(!this->decompress(block.meta_hot_container)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta hot!" << std::endl;
			}
		}

		if(block.meta_cold_container.getSizeCompressed()){
			if(!this->decompress(block.meta_cold_container)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta cold!" << std::endl;
			}
		}

		if(block.gt_rle_container.getSizeCompressed()){
			if(!this->decompress(block.gt_rle_container)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT RLE!" << std::endl;
			}
		}

		if(block.gt_simple_container.getSizeCompressed()){
			if(!this->decompress(block.gt_simple_container)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT other!" << std::endl;
			}
		}

		if(block.gt_support_data_container.getSizeCompressed()){
			if(!this->decompress(block.gt_support_data_container)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT support!" << std::endl;
			}
		}

		if(block.meta_info_map_ids.getSizeCompressed()){
			if(!this->decompress(block.meta_info_map_ids)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to meta INFO ids!" << std::endl;
			}
		}

		if(block.meta_filter_map_ids.getSizeCompressed()){
			if(!this->decompress(block.meta_filter_map_ids)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta FILTER ids!" << std::endl;
			}
		}

		if(block.meta_format_map_ids.getSizeCompressed()){
			if(!this->decompress(block.meta_format_map_ids)){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta FORMAT ids!" << std::endl;
			}
		}

		for(U32 i = 0; i < block.header.n_info_streams; ++i){
			if(block.info_containers[i].getSizeCompressed()){
				if(!this->decompress(block.info_containers[i])){
					std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO " << i << std::endl;
				}
			}
		}

		for(U32 i = 0; i < block.header.n_format_streams; ++i){
			if(block.format_containers[i].getSizeCompressed()){
				if(!this->decompress(block.format_containers[i])){
					std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress FORMAT " << i << std::endl;
				}
			}
		}
		return true;
	}

	bool decompress(algorithm::PermutationManager& permutation_manager){
		return(this->zstd_codec.decompress(permutation_manager));
	}

	/**<
	 * Decompress an abstract data container
	 * @param container Target container
	 * @return          Returns TRUE upon success or FALSE otherwise
	 */
	bool decompress(container_type& container){
		if(container.header.data_header.controller.encoder == YON_ENCODE_ZSTD){
			if(!this->zstd_codec.decompress(container)){ std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress!" << std::endl; return false; }
		} else if(container.header.data_header.controller.encoder == YON_ENCODE_NONE){
			if(!this->no_codec.decompress(container)){ std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress!" << std::endl; return false; }
		} else if(container.header.data_header.controller.encoder == YON_ENCODE_ZPAQ){
			if(!this->zpaq_codec.decompress(container)){ std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "Failed to decompress!" << std::endl; return false; }
		} else {
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
			return false;
		}

		if(container.header.data_header.controller.mixedStride){
			if(container.header.stride_header.controller.encoder == YON_ENCODE_ZSTD){
				this->zstd_codec.decompressStrides(container);
			} else if (container.header.stride_header.controller.encoder == YON_ENCODE_NONE){
				this->no_codec.decompressStrides(container);
			} else if (container.header.stride_header.controller.encoder == YON_ENCODE_ZPAQ){
				this->zpaq_codec.decompressStrides(container);
			} else {
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
				return false;
			}
		}
		return true;
	}

public:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
	zpaq_codec_type zpaq_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
