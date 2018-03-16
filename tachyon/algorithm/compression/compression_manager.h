#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "../../containers/variantblock.h"
#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "zpaq_codec.h"

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
		if(block.meta_contig_container.n_entries)     zstd_codec.compress(block.meta_contig_container);
		if(block.meta_positions_container.n_entries)  zstd_codec.compress(block.meta_positions_container);
		if(block.meta_refalt_container.n_entries)     zstd_codec.compress(block.meta_refalt_container);
		if(block.meta_controller_container.n_entries) zstd_codec.compress(block.meta_controller_container);
		if(block.meta_quality_container.n_entries)    zstd_codec.compress(block.meta_quality_container);
		if(block.meta_names_container.n_entries){
			//zpaq_codec.compress(block.meta_names_container, false);
			zstd_codec.compress(block.meta_names_container);
		}

		const std::string zpaq_cmd = "4";
		if(block.gt_rle8_container.n_entries){
			zstd_codec.compress(block.gt_rle8_container);
			//zstd_codec.compressStrides(block.gt_rle8_container);
		}
		if(block.gt_rle16_container.n_entries){
			zstd_codec.compress(block.gt_rle16_container);
			//zstd_codec.compressStrides(block.gt_rle16_container);
		}
		if(block.gt_rle32_container.n_entries){
			zstd_codec.compress(block.gt_rle32_container);
			//zstd_codec.compressStrides(block.gt_rle32_container);
		}
		if(block.gt_rle64_container.n_entries){
			zstd_codec.compress(block.gt_rle64_container);
			//zstd_codec.compressStrides(block.gt_rle64_container);
		}

		if(block.meta_alleles_container.n_entries){
			zstd_codec.compress(block.meta_alleles_container);
			//zstd_codec.compressStrides(block.meta_alleles_container);
		}

		if(block.gt_simple_container.n_entries)       zstd_codec.compress(block.gt_simple_container);
		if(block.gt_support_data_container.n_entries) zstd_codec.compress(block.gt_support_data_container);
		if(block.meta_info_map_ids.n_entries)         zstd_codec.compress(block.meta_info_map_ids);
		if(block.meta_filter_map_ids.n_entries)       zstd_codec.compress(block.meta_filter_map_ids);
		if(block.meta_format_map_ids.n_entries)       zstd_codec.compress(block.meta_format_map_ids);

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(20);
				//zpaq_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
			}
			else {
				zstd_codec.setCompressionLevel(20);
				//zstd_codec.compress(block.info_containers[i]);
				//zpaq_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
			}
			//zstd_codec.compress(block.info_containers[i]);
			//zpaq_codec.compress(block.info_containers[i], "4");
			zstd_codec.compress(block.info_containers[i]);
		}

		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
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

		if(block.meta_contig_container.getSizeCompressed())     this->decompress(block.meta_contig_container);
		if(block.meta_positions_container.getSizeCompressed())  this->decompress(block.meta_positions_container);
		if(block.meta_refalt_container.getSizeCompressed())     this->decompress(block.meta_refalt_container);
		if(block.meta_alleles_container.getSizeCompressed())    this->decompress(block.meta_alleles_container);
		if(block.meta_controller_container.getSizeCompressed()) this->decompress(block.meta_controller_container);
		if(block.meta_quality_container.getSizeCompressed())    this->decompress(block.meta_quality_container);
		if(block.gt_simple_container.getSizeCompressed())       this->decompress(block.gt_simple_container);
		if(block.gt_support_data_container.getSizeCompressed()) this->decompress(block.gt_support_data_container);
		if(block.meta_info_map_ids.getSizeCompressed())         this->decompress(block.meta_info_map_ids);
		if(block.meta_names_container.getSizeCompressed())      this->decompress(block.meta_names_container);
		if(block.meta_filter_map_ids.getSizeCompressed())       this->decompress(block.meta_filter_map_ids);
		if(block.meta_format_map_ids.getSizeCompressed())       this->decompress(block.meta_format_map_ids);
		if(block.gt_rle8_container.getSizeCompressed())         this->decompress(block.gt_rle8_container);
		if(block.gt_rle16_container.getSizeCompressed())        this->decompress(block.gt_rle16_container);
		if(block.gt_rle32_container.getSizeCompressed())        this->decompress(block.gt_rle32_container);
		if(block.gt_rle64_container.getSizeCompressed())        this->decompress(block.gt_rle64_container);

		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
			if(block.format_containers[i].getSizeCompressed())
				this->decompress(block.format_containers[i]);
		}

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(block.info_containers[i].getSizeCompressed())
				this->decompress(block.info_containers[i]);
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
		//std::cout.write(container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());

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
			//std::cout.write(container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
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
