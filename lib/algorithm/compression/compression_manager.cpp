#include "compression_manager.h"

namespace tachyon{
namespace algorithm{

bool CompressionManager::Compress(variant_block_type& block, const BYTE general_level, const BYTE float_level){
	zstd_codec.SetCompressionLevel(general_level);
	zstd_codec.SetCompressionLevelData(float_level);

	if(block.header.controller.hasGTPermuted)
		zstd_codec.Compress(block.ppa_manager);

	zstd_codec.SetCompressionLevel(general_level);

	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].header.n_entries){
			zstd_codec.Compress(block.base_containers[i]);
			//std::cerr << "Compress: " << i << ": " << block.base_containers[i].buffer_data_uncompressed.size() << "->" << block.base_containers[i].buffer_data.size() << std::endl;
		}
	}

	for(U32 i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
		   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
			zstd_codec.SetCompressionLevelData(float_level);
			zstd_codec.SetCompressionLevelStrides(general_level);
		}
		else {
			zstd_codec.SetCompressionLevel(general_level);
		}
		zstd_codec.Compress(block.info_containers[i]);
		//std::cerr << "Compress INFO: " << i << ": " << block.info_containers[i].buffer_data_uncompressed.size() << "->" << block.info_containers[i].buffer_data.size() << std::endl;
	}

	for(U32 i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
		   block.format_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
			zstd_codec.SetCompressionLevelData(float_level);
			zstd_codec.SetCompressionLevelStrides(general_level);
		}
		else {
			zstd_codec.SetCompressionLevel(general_level);
		}
		zstd_codec.Compress(block.format_containers[i]);
		//std::cerr << "Compress FORMAT: " << i << ": " << block.format_containers[i].buffer_data_uncompressed.size() << "->" << block.format_containers[i].buffer_data.size() << std::endl;
	}

	return true;
}

bool CompressionManager::Decompress(variant_block_type& block){
	if(block.ppa_manager.PPA.size()){
		if(!this->Decompress(block.ppa_manager)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT permutation information!" << std::endl;
			return false;
		}
	}

	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].GetSizeCompressed()){
			if(!this->Decompress(block.base_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress basic container!" << std::endl;
				return false;
			}
		}
	}

	for(U32 i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].GetSizeCompressed()){
			if(!this->Decompress(block.info_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO container " << i << "/" << block.footer.n_info_streams << "!" << std::endl;
				return false;
			}
		}
	}

	for(U32 i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].GetSizeCompressed()){
			if(!this->Decompress(block.format_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress FORMAT container " << i << "/" << block.footer.n_format_streams << "!" << std::endl;
				return false;
			}
		}
	}

	return true;
}

bool CompressionManager::Decompress(algorithm::PermutationManager& permutation_manager){
	return(this->zstd_codec.Decompress(permutation_manager));
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
