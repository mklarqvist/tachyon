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
		const BYTE zstd_compression_level = 19;
		zstd_codec.setCompressionLevel(zstd_compression_level);
		zstd_codec.setCompressionLevelData(3);
		if(block.header.controller.hasGTPermuted)            zstd_codec.compress(block.ppa_manager);
		zstd_codec.setCompressionLevel(zstd_compression_level);
		if(block.meta_contig_container.header.n_entries)     zstd_codec.compress(block.meta_contig_container);
		if(block.meta_positions_container.header.n_entries)  zstd_codec.compress(block.meta_positions_container);
		if(block.meta_refalt_container.header.n_entries)     zstd_codec.compress(block.meta_refalt_container);
		if(block.meta_controller_container.header.n_entries) zstd_codec.compress(block.meta_controller_container);
		if(block.meta_quality_container.header.n_entries)    zstd_codec.compress(block.meta_quality_container);
		if(block.meta_names_container.header.n_entries){
			//zpaq_codec.compress(block.meta_names_container, false);
			zstd_codec.compress(block.meta_names_container);
		}

		//const std::string zpaq_cmd = "4";
		if(block.gt_rle8_container.header.n_entries){
			zstd_codec.compress(block.gt_rle8_container);
			//zstd_codec.compressStrides(block.gt_rle8_container);
		}
		if(block.gt_rle16_container.header.n_entries){
			zstd_codec.compress(block.gt_rle16_container);
			//zstd_codec.compressStrides(block.gt_rle16_container);
		}
		if(block.gt_rle32_container.header.n_entries){
			zstd_codec.compress(block.gt_rle32_container);
			//zstd_codec.compressStrides(block.gt_rle32_container);
		}
		if(block.gt_rle64_container.header.n_entries){
			zstd_codec.compress(block.gt_rle64_container);
			//zstd_codec.compressStrides(block.gt_rle64_container);
		}

		if(block.meta_alleles_container.header.n_entries){
			zstd_codec.compress(block.meta_alleles_container);
			//zstd_codec.compressStrides(block.meta_alleles_container);
		}

		if(block.gt_simple8_container.header.n_entries)      zstd_codec.compress(block.gt_simple8_container);
		if(block.gt_simple16_container.header.n_entries)     zstd_codec.compress(block.gt_simple16_container);
		if(block.gt_simple32_container.header.n_entries)     zstd_codec.compress(block.gt_simple32_container);
		if(block.gt_simple64_container.header.n_entries)     zstd_codec.compress(block.gt_simple64_container);
		if(block.gt_support_data_container.header.n_entries) zstd_codec.compress(block.gt_support_data_container);
		if(block.meta_info_map_ids.header.n_entries)         zstd_codec.compress(block.meta_info_map_ids);
		if(block.meta_filter_map_ids.header.n_entries)       zstd_codec.compress(block.meta_filter_map_ids);
		if(block.meta_format_map_ids.header.n_entries)       zstd_codec.compress(block.meta_format_map_ids);

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.setCompressionLevelData(3);
				zstd_codec.setCompressionLevelStrides(zstd_compression_level);
				//zpaq_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
				//zstd_codec.compress(block.info_containers[i]);
			}
			else {
				zstd_codec.setCompressionLevel(zstd_compression_level);
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
				zstd_codec.setCompressionLevelStrides(zstd_compression_level);
				//zpaq_codec.compress(block.format_containers[i]);
				//zstd_codec.compress(block.format_containers[i]);
			}
			else {
				zstd_codec.setCompressionLevel(zstd_compression_level);
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
				return false;
			}
		}

		if(block.meta_contig_container.getSizeCompressed())     if(!this->decompress(block.meta_contig_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta contig information!" << std::endl; return false; }
		if(block.meta_positions_container.getSizeCompressed())  if(!this->decompress(block.meta_positions_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta positions information!" << std::endl; return false; }
		if(block.meta_refalt_container.getSizeCompressed())     if(!this->decompress(block.meta_refalt_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta ref_alt information!" << std::endl; return false; }
		if(block.meta_alleles_container.getSizeCompressed())    if(!this->decompress(block.meta_alleles_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta alleles information!" << std::endl; return false; }
		if(block.meta_controller_container.getSizeCompressed()) if(!this->decompress(block.meta_controller_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta controller information!" << std::endl; return false; }
		if(block.meta_quality_container.getSizeCompressed())    if(!this->decompress(block.meta_quality_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta quality information!" << std::endl; return false; }
		if(block.gt_support_data_container.getSizeCompressed()) if(!this->decompress(block.gt_support_data_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotype support information!" << std::endl; return false; }
		if(block.meta_info_map_ids.getSizeCompressed())         if(!this->decompress(block.meta_info_map_ids)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta INFO maps information!" << std::endl; return false; }
		if(block.meta_names_container.getSizeCompressed())      if(!this->decompress(block.meta_names_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta names information!" << std::endl; return false; }
		if(block.meta_filter_map_ids.getSizeCompressed())       if(!this->decompress(block.meta_filter_map_ids)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta FILTER maps information!" << std::endl; return false; }
		if(block.meta_format_map_ids.getSizeCompressed())       if(!this->decompress(block.meta_format_map_ids)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress meta FORMAT maps information!" << std::endl; return false; }
		if(block.gt_rle8_container.getSizeCompressed())         if(!this->decompress(block.gt_rle8_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (RLE-8) information!" << std::endl; return false; }
		if(block.gt_rle16_container.getSizeCompressed())        if(!this->decompress(block.gt_rle16_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (RLE-16) information!" << std::endl; return false; }
		if(block.gt_rle32_container.getSizeCompressed())        if(!this->decompress(block.gt_rle32_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (RLE-32) information!" << std::endl; return false; }
		if(block.gt_rle64_container.getSizeCompressed())        if(!this->decompress(block.gt_rle64_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (RLE-64) information!" << std::endl; return false; }
		if(block.gt_simple8_container.getSizeCompressed())      if(!this->decompress(block.gt_simple8_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (simple-8) information!" << std::endl; return false; }
		if(block.gt_simple16_container.getSizeCompressed())     if(!this->decompress(block.gt_simple16_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (simple-16) information!" << std::endl; return false; }
		if(block.gt_simple32_container.getSizeCompressed())     if(!this->decompress(block.gt_simple32_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (simple-32) information!" << std::endl; return false; }
		if(block.gt_simple64_container.getSizeCompressed())     if(!this->decompress(block.gt_simple64_container)){ std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress genotypes (simple-64) information!" << std::endl; return false; }

		for(U32 i = 0; i < block.footer.n_info_streams; ++i){
			if(block.info_containers[i].getSizeCompressed()){
				if(!this->decompress(block.info_containers[i])){
					std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO container " << i << "/" << block.footer.n_info_streams << "!" << std::endl;
					return false;
				}
			}
		}

		for(U32 i = 0; i < block.footer.n_format_streams; ++i){
			if(block.format_containers[i].getSizeCompressed()){
				if(!this->decompress(block.format_containers[i])){
					std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress FORMAT container " << i << "/" << block.footer.n_format_streams << "!" << std::endl;
					return false;
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
			if(!this->zstd_codec.decompress(container)){
				std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress data!" << std::endl;
				return false;
			}
		} else if(container.header.data_header.controller.encoder == YON_ENCODE_NONE){
			if(!this->no_codec.decompress(container)){
				std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress data!" << std::endl;
				return false;
			}
		} else if(container.header.data_header.controller.encoder == YON_ENCODE_ZPAQ){
			if(!this->zpaq_codec.decompress(container)){
				std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "Failed to decompress data!" << std::endl;
				return false;
			}
		} else {
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
			return false;
		}

		if(container.header.data_header.controller.mixedStride){
			if(container.header.stride_header.controller.encoder == YON_ENCODE_ZSTD){
				if(!this->zstd_codec.decompressStrides(container)){ std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress strides!" << std::endl; return false; }
			} else if (container.header.stride_header.controller.encoder == YON_ENCODE_NONE){
				if(!this->no_codec.decompressStrides(container)){ std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress strides!" << std::endl; return false; }
			} else if (container.header.stride_header.controller.encoder == YON_ENCODE_ZPAQ){
				if(!this->zpaq_codec.decompressStrides(container)){ std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "Failed to decompress strides!" << std::endl; return false; }
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
