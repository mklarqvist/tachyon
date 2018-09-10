#include "core/genotypes.h"
#include "compression_manager.h"

#include "algorithm/compression/fastdelta.h"

#include "containers/stride_container.h"

namespace tachyon {
namespace algorithm {

bool CompressionManager::Compress(variant_block_type& block,
                                  const uint8_t general_level)
{
	zstd_codec.SetCompressionLevel(general_level);
	zstd_codec.SetCompressionLevelData(general_level);

	if(block.header.controller.has_gt_permuted){
		zstd_codec.SetCompressionLevel(22);
		zstd_codec.Compress(block.base_containers[YON_BLK_PPA], *block.gt_ppa);
		//this->CompressEvaluate(block.base_containers[YON_BLK_PPA], block.gt_ppa);
	}

	zstd_codec.SetCompressionLevel(general_level);

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		/*
		if(block.base_containers[i].data_uncompressed.size()){
			if(block.base_containers[i].header.data_header.IsUniform() == false &&
			   block.base_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				if(this->memory_basic[i].permanent_codec >= 0){
					//std::cerr << "using permanent base" << std::endl;
					switch(this->memory_basic[i].codec){
						case(0): zstd_codec.Compress(block.base_containers[i]); break;
						case(1): this->CompressCodec1(block.base_containers[i]); break;
						case(2): this->CompressCodec2(block.base_containers[i]); break;
						case(3): this->CompressCodec3(block.base_containers[i]); break;
						case(4): this->CompressCodec4(block.base_containers[i]); break;
						case(5): this->CompressDelta(block.base_containers[i]); break;
					}
				} else if(this->memory_basic[i].times_seen == 10){
					//std::cerr << utility::timestamp("DEBUG") << "eval base " << YON_BLK_PRINT_NAMES[i] << std::endl;
					this->CompressEvaluate(block.base_containers[i], this->memory_basic, i);
					this->memory_basic[i].times_seen = 0;
				} else {
					//std::cerr << "reuse: " << (int)this->memory_basic[i].codec << std::endl;
					switch(this->memory_basic[i].codec){
						case(0): zstd_codec.Compress(block.base_containers[i]); break;
						case(1): this->CompressCodec1(block.base_containers[i]); break;
						case(2): this->CompressCodec2(block.base_containers[i]); break;
						case(3): this->CompressCodec3(block.base_containers[i]); break;
						case(4): this->CompressCodec4(block.base_containers[i]); break;
						case(5): this->CompressDelta(block.base_containers[i]); break;
					}
					++this->memory_basic[i].times_seen;
				}
			}  else {
				zstd_codec.Compress(block.base_containers[i]);
			}
		}
		*/
			zstd_codec.Compress(block.base_containers[i]);
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
		   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
			zstd_codec.SetCompressionLevelData(general_level);
			zstd_codec.SetCompressionLevelStrides(general_level);
		}
		else {
			zstd_codec.SetCompressionLevel(general_level);
		}

		if(block.info_containers[i].header.n_entries){
			/*
			if(block.info_containers[i].header.data_header.IsUniform() == false &&
			   block.info_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				const int32_t global_id = block.footer.info_offsets[i].GetGlobalKey();
				//std::cerr << "info global: " << global_id << std::endl;
				if(this->memory_info[global_id].permanent_codec >= 0){
					//std::cerr << "using permanent info" << std::endl;
					switch(this->memory_info[global_id].codec){
						case(0): zstd_codec.Compress(block.info_containers[i]); break;
						case(1): this->CompressCodec1(block.info_containers[i]); break;
						case(2): this->CompressCodec2(block.info_containers[i]); break;
						case(3): this->CompressCodec3(block.info_containers[i]); break;
						case(4): this->CompressCodec4(block.info_containers[i]); break;
						case(5): this->CompressDelta(block.info_containers[i]); break;
					}
				} else if(this->memory_info[global_id].times_seen == 10){
					//std::cerr << utility::timestamp("DEBUG") << "eval info " << global_id << std::endl;
					this->CompressEvaluate(block.info_containers[i], this->memory_info, global_id);
					this->memory_info[global_id].times_seen = 0;
				} else {
					//std::cerr << "reuse: " << (int)this->memory_info[i].codec << std::endl;
					switch(this->memory_info[global_id].codec){
						case(0): zstd_codec.Compress(block.info_containers[i]); break;
						case(1): this->CompressCodec1(block.info_containers[i]); break;
						case(2): this->CompressCodec2(block.info_containers[i]); break;
						case(3): this->CompressCodec3(block.info_containers[i]); break;
						case(4): this->CompressCodec4(block.info_containers[i]); break;
						case(5): this->CompressDelta(block.info_containers[i]); break;
					}
					++this->memory_info[global_id].times_seen;
				}
			} else {
				zstd_codec.Compress(block.info_containers[i]);
			}
			*/
			std::cerr << "Info mixed stride: " << block.info_containers[i].header.data_header.HasMixedStride() << std::endl;
			zstd_codec.Compress(block.info_containers[i]);
		}
	}

	//io::BasicBuffer temp;
	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].header.n_entries){
			if(block.format_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.format_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.SetCompressionLevelData(general_level);
				zstd_codec.SetCompressionLevelStrides(general_level);
			}
			else {
				zstd_codec.SetCompressionLevel(general_level);
			}

			if(block.format_containers[i].header.n_entries){
				/*
				const int32_t global_id = block.footer.format_offsets[i].GetGlobalKey();
				if(block.format_containers[i].header.data_header.IsUniform() == false &&
				   block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
				{
					if(this->memory_format[global_id].permanent_codec >= 0){
						//std::cerr << "using permanent format" << std::endl;
						switch(this->memory_format[global_id].codec){
							case(0): zstd_codec.Compress(block.format_containers[i]); break;
							case(1): this->CompressCodec1(block.format_containers[i]); break;
							case(2): this->CompressCodec2(block.format_containers[i]); break;
							case(3): this->CompressCodec3(block.format_containers[i]); break;
							case(4): this->CompressCodec4(block.format_containers[i]); break;
							case(5): this->CompressDelta(block.format_containers[i]); break;
						}
					} else if(this->memory_format[global_id].times_seen == 10){
						//std::cerr << utility::timestamp("DEBUG") << "eval format " << global_id << std::endl;
						this->CompressEvaluate(block.format_containers[i], this->memory_format, global_id);
						this->memory_format[global_id].times_seen = 0;
					} else {
						//std::cerr << "reuse: " << (int)this->memory_format[i].codec << std::endl;
						switch(this->memory_format[global_id].codec){
							case(0): zstd_codec.Compress(block.format_containers[i]); break;
							case(1): this->CompressCodec1(block.format_containers[i]); break;
							case(2): this->CompressCodec2(block.format_containers[i]); break;
							case(3): this->CompressCodec3(block.format_containers[i]); break;
							case(4): this->CompressCodec4(block.format_containers[i]); break;
							case(5): this->CompressDelta(block.format_containers[i]); break;
						}
						++this->memory_format[global_id].times_seen;
					}
				} else {
					zstd_codec.Compress(block.format_containers[i]);
				}
				*/
				/*
				if(block.format_containers[i].header.data_header.IsUniform() == false &&
				   (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B || block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_16B || block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_8B))
				{
					temp.resize(block.format_containers[i].data_uncompressed.size() + 65536);
					containers::StrideContainer s(block.format_containers[i]);

					uint8_t byte_size =

					for(int k = 0; k < s.size(); ++s){
						block.format_containers[i].data_uncompressed
					}
				}
				*/
				zstd_codec.Compress(block.format_containers[i]);
			}
		}
	}
	return true;
}

bool CompressionManager::Decompress(variant_block_type& block){
	if(block.base_containers[YON_BLK_PPA].GetSizeCompressed() &&
	   block.base_containers[YON_BLK_PPA].header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
	{
		if(!this->Decompress(block.base_containers[YON_BLK_PPA], *block.gt_ppa)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT permutation information!" << std::endl;
			return false;
		}
	}

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].GetSizeCompressed() &&
		   block.base_containers[i].IsEncrypted() == false)
		{
			if(!this->Decompress(block.base_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress basic container!" << std::endl;
				return false;
			}
		}
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].GetSizeCompressed() &&
		   block.info_containers[i].IsEncrypted() == false)
		{
			if(!this->Decompress(block.info_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO container " << i << "/" << block.footer.n_info_streams << "!" << std::endl;
				return false;
			}
		}
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].GetSizeCompressed() &&
		   block.format_containers[i].IsEncrypted() == false)
		{
			if(!this->Decompress(block.format_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress FORMAT container " << i << "/" << block.footer.n_format_streams << "!" << std::endl;
				return false;
			}
		}
	}

	return true;
}

bool CompressionManager::Decompress(container_type& container, yon_gt_ppa& gt_ppa){
	return(this->zstd_codec.Decompress(container, gt_ppa));
}

bool CompressionManager::Decompress(container_type& container){
	// Ascertain that data is not encrypted
	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if(container.header.data_header.controller.encoder == YON_ENCODE_ZSTD){
		if(!this->zstd_codec.Decompress(container)){
			std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress data!" << std::endl;
			return false;
		}
	} else if(container.header.data_header.controller.encoder == YON_ENCODE_NONE){
		if(!this->no_codec.Decompress(container)){
			std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress data!" << std::endl;
			return false;
		}
	} else if(container.header.data_header.controller.encoder == YON_ENCODE_ZPAQ){
		std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "ZPAQ is no longer supported!" << std::endl;
		return false;
	} else {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
		return false;
	}

	if(container.header.data_header.controller.mixedStride){
		if(container.header.stride_header.controller.encoder == YON_ENCODE_ZSTD){
			if(!this->zstd_codec.DecompressStrides(container)){ std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress strides!" << std::endl; return false; }
		} else if (container.header.stride_header.controller.encoder == YON_ENCODE_NONE){
			if(!this->no_codec.DecompressStrides(container)){ std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress strides!" << std::endl; return false; }
		} else if (container.header.stride_header.controller.encoder == YON_ENCODE_ZPAQ){
			std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "ZPAQ is no longer supported!" << std::endl;
			return false;
		} else {
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
			return false;
		}
		//std::cout.write(container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
	}
	return true;
}

}
}
