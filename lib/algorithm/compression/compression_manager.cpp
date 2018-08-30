#include "core/genotypes.h"
#include "compression_manager.h"

#include "codecfactory.h"
#include "deltautil.h"

namespace tachyon{
namespace algorithm{

bool CompressionManager::Compress(variant_block_type& block, const uint8_t general_level, const uint8_t float_level){
	zstd_codec.SetCompressionLevel(general_level);
	zstd_codec.SetCompressionLevelData(float_level);

	if(block.header.controller.hasGTPermuted){
		zstd_codec.SetCompressionLevel(22);
		zstd_codec.Compress(block.base_containers[YON_BLK_PPA], *block.gt_ppa);
	}

	zstd_codec.SetCompressionLevel(general_level);

	FastPForLib::IntegerCODEC* codec1 = new FastPForLib::CompositeCodec<FastPForLib::SIMDOPTPFor<4, FastPForLib::Simple16<false>>, FastPForLib::VariableByte>();
	FastPForLib::IntegerCODEC* codec2 = new MaskedVByte();
	FastPForLib::IntegerCODEC* codec3 = new FastPForLib::StreamVByte();
	FastPForLib::IntegerCODEC* codec4 = new FastPForLib::VariableByte();

	io::BasicBuffer temp;
	io::BasicBuffer backup;
	io::BasicBuffer c1,c2,c3,c4;

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].header.n_entries){
			if(block.base_containers[i].header.data_header.IsUniform() == false &&
			   block.base_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				temp.reset();
				assert(block.base_containers[i].data_uncompressed.size() % sizeof(int32_t) == 0);
				temp.resize(block.base_containers[i].data_uncompressed.size() + 1000*sizeof(int32_t));
				size_t compressedsize = temp.capacity() / sizeof(int32_t);

				// Move uncompressed data into backup
				backup = std::move(block.base_containers[i].data_uncompressed);

				const uint32_t* data = reinterpret_cast<const uint32_t*>(backup.data());
				uint32_t* dst_data   = reinterpret_cast<uint32_t*>(temp.data());
				assert(block.base_containers[i].data_uncompressed.size() % sizeof(int32_t) == 0);

				// Compress data into temp buffer
				codec1->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.base_containers[i].data_uncompressed = temp;
				// Compress
				zstd_codec.Compress(block.base_containers[i]);
				// Copy compressed data
				c1 = std::move(block.base_containers[i].data);
				block.base_containers[i].data.resize(65536);
				// Print
				std::cerr << "preprocessor " << i << " simdoptpfor: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c1.size() << std::endl;


				// Codec2
				temp.reset();
				compressedsize = temp.capacity()/sizeof(int32_t);

				// Compress data into temp buffer
				codec2->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.base_containers[i].data_uncompressed = temp;
				// Copy compressed data
				zstd_codec.Compress(block.base_containers[i]);
				// Copy compressed data
				c2 = std::move(block.base_containers[i].data);
				block.base_containers[i].data.resize(65536);
				// Print
				std::cerr << "preprocessor " << i << " maskedvbyte: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c2.size() << std::endl;

				// Codec3
				temp.reset();
				compressedsize = temp.capacity()/sizeof(int32_t);

				// Compress data into temp buffer
				codec3->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.base_containers[i].data_uncompressed = temp;
				// Copy compressed data
				zstd_codec.Compress(block.base_containers[i]);
				// Copy compressed data
				c3 = std::move(block.base_containers[i].data);
				block.base_containers[i].data.resize(65536);
				// Print
				std::cerr << "preprocessor " << i << " streamvbyte: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c3.size() << std::endl;

				// Move back
				block.base_containers[i].data_uncompressed = std::move(backup);
			}

			zstd_codec.Compress(block.base_containers[i]);

			if(block.base_containers[i].header.data_header.IsUniform() == false &&
			   block.base_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				std::cerr << "Base Compress: " << i << ": " << block.base_containers[i].data_uncompressed.size() << "->" << block.base_containers[i].data.size() << " for type " << block.base_containers[i].header.data_header.GetPrimitiveType() << std::endl;

				uint32_t codec_min = block.base_containers[i].data.size();
				uint8_t codec_chosen = 0;
				if(c1.size() < codec_min && c1.size() != 0){ codec_min = c1.size(); codec_chosen = 1; }
				if(c2.size() < codec_min && c2.size() != 0){ codec_min = c2.size(); codec_chosen = 2; }
				if(c3.size() < codec_min && c3.size() != 0){ codec_min = c3.size(); codec_chosen = 3; }

				std::cerr << "Base Winner: " << (int)codec_chosen << " with " << codec_min << " difference " << (float)block.base_containers[i].data.size()/codec_min << std::endl;
				if(codec_chosen == 1)     { block.base_containers[i].header.data_header.cLength = c1.size(); block.base_containers[i].data = std::move(c1); }
				else if(codec_chosen == 2){ block.base_containers[i].header.data_header.cLength = c2.size(); block.base_containers[i].data = std::move(c2); }
				else if(codec_chosen == 3){ block.base_containers[i].header.data_header.cLength = c3.size(); block.base_containers[i].data = std::move(c3); }
			}

		}
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
		   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
			zstd_codec.SetCompressionLevelData(float_level);
			zstd_codec.SetCompressionLevelStrides(general_level);
		}
		else {
			zstd_codec.SetCompressionLevel(general_level);
		}
		block.info_containers[i].ReformatInteger();
		zstd_codec.Compress(block.info_containers[i]);
		//std::cerr << "Compress INFO: " << i << ": " << block.info_containers[i].buffer_data_uncompressed.size() << "->" << block.info_containers[i].buffer_data.size() << std::endl;
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].header.n_entries){
			if(block.format_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.format_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE){
				zstd_codec.SetCompressionLevelData(float_level);
				zstd_codec.SetCompressionLevelStrides(general_level);
			}
			else {
				zstd_codec.SetCompressionLevel(general_level);
			}


			if(block.format_containers[i].header.data_header.IsUniform() == false &&
			   block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				temp.reset();
				temp.resize(block.format_containers[i].data_uncompressed.size() + 1000*sizeof(int32_t));
				size_t compressedsize = temp.capacity()/sizeof(int32_t);

				// Move uncompressed data into backup
				backup = std::move(block.format_containers[i].data_uncompressed);

				const uint32_t* data = reinterpret_cast<const uint32_t*>(backup.data());
				uint32_t* dst_data   = reinterpret_cast<uint32_t*>(temp.data());
				assert(block.format_containers[i].data_uncompressed.size() % sizeof(int32_t) == 0);

				// Compress data into temp buffer
				codec1->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp

				block.format_containers[i].data_uncompressed = temp;
				// Compress
				zstd_codec.Compress(block.format_containers[i]);
				// Copy compressed data
				c1 = std::move(block.format_containers[i].data);
				block.format_containers[i].data.reset();
				// Print
				std::cerr << "preprocessor " << i << " simdoptpfor: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c1.size() << std::endl;

				// Codec2
				temp.reset();
				compressedsize = temp.capacity()/sizeof(int32_t);

				// Compress data into temp buffer
				codec2->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.format_containers[i].data_uncompressed = temp;
				// Copy compressed data
				zstd_codec.Compress(block.format_containers[i]);
				// Copy compressed data
				c2 = std::move(block.format_containers[i].data);
				block.format_containers[i].data.reset();
				// Print
				std::cerr << "preprocessor " << i << " maskedvbyte: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c2.size() << std::endl;

				// Codec3
				temp.reset();
				compressedsize = temp.capacity()/sizeof(int32_t);

				// Compress data into temp buffer
				codec3->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.format_containers[i].data_uncompressed = temp;
				// Copy compressed data
				zstd_codec.Compress(block.format_containers[i]);
				// Copy compressed data
				c3 = std::move(block.format_containers[i].data);
				block.format_containers[i].data.reset();
				// Print
				std::cerr << "preprocessor " << i << " streamvbyte: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c3.size() << std::endl;

				// Codec4
				temp.reset();
				compressedsize = temp.capacity()/sizeof(int32_t);


				// Compress data into temp buffer
				codec4->encodeArray(
						data,
						backup.size()/sizeof(int32_t),
						dst_data,
						compressedsize);

				temp.n_chars_ = compressedsize*sizeof(int32_t);
				// Set uncompressed data to temp
				block.format_containers[i].data_uncompressed = temp;
				// Copy compressed data
				zstd_codec.Compress(block.format_containers[i]);
				// Copy compressed data
				c4 = std::move(block.format_containers[i].data);
				block.format_containers[i].data.reset();
				// Print
				std::cerr << "preprocessor " << i << " varint: " << backup.size() << "->" << compressedsize*sizeof(int32_t) << "->" << c4.size() << std::endl;

				// Move back
				block.format_containers[i].data_uncompressed = std::move(backup);
			}


			zstd_codec.Compress(block.format_containers[i]);


			if(block.format_containers[i].header.data_header.IsUniform() == false &&
			   block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				std::cerr << "Format Compress: " << i << ": " << block.format_containers[i].data_uncompressed.size() << "->" << block.format_containers[i].data.size() << " for type " << block.format_containers[i].header.data_header.GetPrimitiveType() << std::endl;

				uint32_t codec_min = block.format_containers[i].data.size();
				uint8_t codec_chosen = 0;
				if(c1.size() < codec_min && c1.size() != 0){ codec_min = c1.size(); codec_chosen = 1; }
				if(c2.size() < codec_min && c2.size() != 0){ codec_min = c2.size(); codec_chosen = 2; }
				if(c3.size() < codec_min && c3.size() != 0){ codec_min = c3.size(); codec_chosen = 3; }
				if(c4.size() < codec_min && c4.size() != 0){ codec_min = c4.size(); codec_chosen = 4; }

				std::cerr << "Format Winner: " << (int)codec_chosen << " with " << codec_min << " difference " << (float)block.format_containers[i].data.size()/codec_min << std::endl;
				if(codec_chosen == 1)     { block.format_containers[i].header.data_header.cLength = c1.size(); block.format_containers[i].data = std::move(c1); }
				else if(codec_chosen == 2){ block.format_containers[i].header.data_header.cLength = c2.size(); block.format_containers[i].data = std::move(c2); }
				else if(codec_chosen == 3){ block.format_containers[i].header.data_header.cLength = c3.size(); block.format_containers[i].data = std::move(c3); }
				else if(codec_chosen == 4){ block.format_containers[i].header.data_header.cLength = c4.size(); block.format_containers[i].data = std::move(c4); }
			}

		}
	}

	delete codec1, codec2, codec3, codec4;

	return true;
}

bool CompressionManager::Decompress(variant_block_type& block){
	if(block.base_containers[YON_BLK_PPA].GetSizeCompressed()){
		if(!this->Decompress(block.base_containers[YON_BLK_PPA], *block.gt_ppa)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT permutation information!" << std::endl;
			return false;
		}
	}

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].GetSizeCompressed()){
			if(!this->Decompress(block.base_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress basic container!" << std::endl;
				return false;
			}
		}
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].GetSizeCompressed()){
			if(!this->Decompress(block.info_containers[i])){
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO container " << i << "/" << block.footer.n_info_streams << "!" << std::endl;
				return false;
			}
		}
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].GetSizeCompressed()){
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
