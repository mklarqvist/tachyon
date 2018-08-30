#include "core/genotypes.h"
#include "compression_manager.h"

#include "codecfactory.h"
#include "deltautil.h"
#include "algorithm/compression/fastdelta.h"

namespace tachyon{
namespace algorithm{

CompressionManager::CompressionManager() :
	codec1(new FastPForLib::CompositeCodec<FastPForLib::SIMDOPTPFor<4, FastPForLib::Simple16<false>>, FastPForLib::VariableByte>()),
	codec2(new MaskedVByte()),
	codec3(new FastPForLib::StreamVByte()),
	codec4(new FastPForLib::CompositeCodec<FastPForLib::SIMDPFor, FastPForLib::VariableByte>())
{

}
CompressionManager::~CompressionManager(){
	delete codec1, codec2, codec3, codec4;
}

int32_t CompressionManager::CompressCodecWrapper(FastPForLib::IntegerCODEC* codec, container_type& container, const bool move){
	if(container.header.data_header.IsUniform() == true ||
	   container.header.data_header.GetPrimitiveType() != YON_TYPE_32B)
	{
		return 0;
	}

	this->support_buffer.reset();
	this->support_buffer.resize(container.data_uncompressed.size() + 65536);
	const uint32_t* data = reinterpret_cast<const uint32_t*>(container.data_uncompressed.data());
	uint32_t* dst_data   = reinterpret_cast<uint32_t*>(this->support_buffer.data());

	// Compress data into temp buffer
	const uint32_t n_entries = container.data_uncompressed.size() / sizeof(uint32_t);
	size_t compressed_size = this->support_buffer.capacity() / sizeof(uint32_t);
	codec->encodeArray(
			data,
			n_entries,
			dst_data,
			compressed_size);

	this->support_buffer.n_chars_ = compressed_size;
	this->backup_buffer           = std::move(container.data_uncompressed);
	container.data_uncompressed   = std::move(this->support_buffer);
	zstd_codec.Compress(container);
	container.data_uncompressed = this->backup_buffer;

	return container.data.size();
}

int32_t CompressionManager::CompressCodec1(container_type& container, const bool move){
	return(this->CompressCodecWrapper(this->codec1, container, move));
}

int32_t CompressionManager::CompressCodec2(container_type& container, const bool move){
	return(this->CompressCodecWrapper(this->codec2, container, move));
}
int32_t CompressionManager::CompressCodec3(container_type& container, const bool move){
	return(this->CompressCodecWrapper(this->codec3, container, move));
}
int32_t CompressionManager::CompressCodec4(container_type& container, const bool move){
	return(this->CompressCodecWrapper(this->codec4, container, move));
}
int32_t CompressionManager::CompressDelta(container_type& container, const bool move){
	if(container.header.data_header.IsUniform() == true ||
	   container.header.data_header.GetPrimitiveType() != YON_TYPE_32B)
	{
		return 0;
	}

	this->support_buffer.reset();
	this->support_buffer.resize(container.data_uncompressed.size() + 65536);
	const uint32_t* data = reinterpret_cast<const uint32_t*>(container.data_uncompressed.data());
	uint32_t* dst_data   = reinterpret_cast<uint32_t*>(this->support_buffer.data());

	// Compress data into temp buffer
	const uint32_t n_entries = container.data_uncompressed.size() / sizeof(uint32_t);

	// Compute deltas
	compute_deltas(data, n_entries, dst_data, data[0]);

	this->support_buffer.n_chars_ = n_entries * sizeof(uint32_t);
	this->backup_buffer           = std::move(container.data_uncompressed);
	container.data_uncompressed   = std::move(this->support_buffer);
	zstd_codec.Compress(container);
	container.data_uncompressed = this->backup_buffer;

	return container.data.size();
}

bool CompressionManager::CompressEvaluate(container_type& container, yon_comp_memory_pair* memory, const uint32_t global_id){
	if(container.header.n_entries == 0) return false;

	uint32_t c_sz[6]; memset(c_sz, 0, sizeof(uint32_t)*6);

	if(container.header.data_header.IsUniform() == true ||
	   container.header.data_header.GetPrimitiveType() != YON_TYPE_32B)
	{
		return false;
	}

	c_sz[1] = this->CompressCodec1(container, false);
	c_sz[2] = this->CompressCodec2(container, false);
	c_sz[3] = this->CompressCodec3(container, false);
	c_sz[4] = this->CompressCodec4(container, false);
	c_sz[5] = this->CompressDelta(container, false);
	zstd_codec.Compress(container);
	c_sz[0] = container.data.size();


	uint32_t codec_min = c_sz[0];
	uint8_t codec_chosen = 0;
	if(c_sz[1] < codec_min && c_sz[1] != 0){ codec_min = c_sz[1]; codec_chosen = 1; }
	if(c_sz[2] < codec_min && c_sz[2] != 0){ codec_min = c_sz[2]; codec_chosen = 2; }
	if(c_sz[3] < codec_min && c_sz[3] != 0){ codec_min = c_sz[3]; codec_chosen = 3; }
	if(c_sz[4] < codec_min && c_sz[4] != 0){ codec_min = c_sz[4]; codec_chosen = 4; }
	if(c_sz[5] < codec_min && c_sz[5] != 0){ codec_min = c_sz[5]; codec_chosen = 5; }

	std::cerr << utility::timestamp("DEBUG") << c_sz[0] << "\t" << c_sz[1] << "\t" << c_sz[2] << "\t" << c_sz[3] << "\t" <<  c_sz[4] << "\t" << c_sz[5] << std::endl;
	std::cerr << "Winner @ " << global_id << ": " << (int)codec_chosen << " with " << codec_min << " difference " << (float)container.data.size()/codec_min << std::endl;
	if(codec_chosen == 1)     {
		container.header.data_header.cLength = c_sz[1];
		this->CompressCodec1(container, true);
		memory[global_id].codec = 1;
	}
	else if(codec_chosen == 2){
		container.header.data_header.cLength = c_sz[2];
		this->CompressCodec2(container, true);
		memory[global_id].codec = 2;
	}
	else if(codec_chosen == 3){
		container.header.data_header.cLength = c_sz[3];
		this->CompressCodec3(container, true);
		memory[global_id].codec = 3;
	}
	else if(codec_chosen == 4){
		container.header.data_header.cLength = c_sz[4];
		this->CompressCodec4(container, true);
		memory[global_id].codec = 4;
	}
	else if(codec_chosen == 5){
		container.header.data_header.cLength = c_sz[5];
		this->CompressDelta(container, true);
		memory[global_id].codec = 5;
	}
	else {
		memory[global_id].codec = 0;
	}

	return true;
}


bool CompressionManager::Compress(variant_block_type& block, const uint8_t general_level, const uint8_t float_level){
	zstd_codec.SetCompressionLevel(general_level);
	zstd_codec.SetCompressionLevelData(float_level);

	if(block.header.controller.hasGTPermuted){
		zstd_codec.SetCompressionLevel(22);
		zstd_codec.Compress(block.base_containers[YON_BLK_PPA], *block.gt_ppa);
	}

	zstd_codec.SetCompressionLevel(general_level);

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].header.n_entries){
			if(block.base_containers[i].header.data_header.IsUniform() == false &&
			   block.base_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				if(this->memory_basic[i].times_seen == 5){
					//std::cerr << "eval" << std::endl;
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
			} else {
				zstd_codec.Compress(block.base_containers[i]);
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

		if(block.info_containers[i].header.n_entries){
			if(block.info_containers[i].header.data_header.IsUniform() == false &&
			   block.info_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				if(this->memory_info[i].times_seen == 5){
					//std::cerr << "eval" << std::endl;
					this->CompressEvaluate(block.info_containers[i], this->memory_info, i);
					this->memory_info[i].times_seen = 0;
				} else {
					//std::cerr << "reuse: " << (int)this->memory_info[i].codec << std::endl;
					switch(this->memory_info[i].codec){
						case(0): zstd_codec.Compress(block.info_containers[i]); break;
						case(1): this->CompressCodec1(block.info_containers[i]); break;
						case(2): this->CompressCodec2(block.info_containers[i]); break;
						case(3): this->CompressCodec3(block.info_containers[i]); break;
						case(4): this->CompressCodec4(block.info_containers[i]); break;
						case(5): this->CompressDelta(block.info_containers[i]); break;
					}
					++this->memory_info[i].times_seen;
				}
			} else {
				zstd_codec.Compress(block.info_containers[i]);
			}
		}
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

			if(block.format_containers[i].header.n_entries){
				if(block.format_containers[i].header.data_header.IsUniform() == false &&
				   block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
				{
					if(this->memory_format[i].times_seen == 5){
						//std::cerr << "eval" << std::endl;
						this->CompressEvaluate(block.format_containers[i], this->memory_format, i);
						this->memory_format[i].times_seen = 0;
					} else {
						//std::cerr << "reuse: " << (int)this->memory_format[i].codec << std::endl;
						switch(this->memory_format[i].codec){
							case(0): zstd_codec.Compress(block.format_containers[i]); break;
							case(1): this->CompressCodec1(block.format_containers[i]); break;
							case(2): this->CompressCodec2(block.format_containers[i]); break;
							case(3): this->CompressCodec3(block.format_containers[i]); break;
							case(4): this->CompressCodec4(block.format_containers[i]); break;
							case(5): this->CompressDelta(block.format_containers[i]); break;
						}
						++this->memory_format[i].times_seen;
					}
				} else {
					zstd_codec.Compress(block.format_containers[i]);
				}
			}
		}
	}
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
