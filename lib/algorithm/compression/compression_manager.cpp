#include "genotypes.h"
#include "compression_manager.h"

namespace tachyon {
namespace algorithm {

bool CompressionManager::Compress(variant_block_type& block,
                                  const uint8_t general_level,
                                  const uint32_t n_samples)
{
	zstd_codec.SetCompressionLevel(general_level);
	zstd_codec.SetCompressionLevelData(general_level);

	if (block.header.controller.has_gt_permuted) {
		zstd_codec.SetCompressionLevel(22);
		zstd_codec.Compress(block.base_containers[YON_BLK_PPA], *block.gt_ppa);
		//this->CompressEvaluate(block.base_containers[YON_BLK_PPA], block.gt_ppa);
	}

	zstd_codec.SetCompressionLevel(general_level);

	for (uint32_t i = 1; i < YON_BLK_N_STATIC; ++i) {
		zstd_codec.Compress(block.base_containers[i]);
	}

	for (uint32_t i = 0; i < block.footer.n_info_streams; ++i) {
		if (block.info_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
		   block.info_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE) {
			zstd_codec.SetCompressionLevelData(general_level);
			zstd_codec.SetCompressionLevelStrides(general_level);
		}
		else {
			zstd_codec.SetCompressionLevel(general_level);
		}

		//std::cerr << "in info-" << i << "/" << block.footer.n_info_streams << " -> " << block.info_containers[i].GetIdx() << "," << block.info_containers[i].GetSizeUncompressed() << ", n=" << block.info_containers[i].header.n_entries << ",ns=" << block.info_containers[i].header.stride_header.stride << std::endl;
		if (block.info_containers[i].header.n_entries) {
			/*
			if (block.info_containers[i].header.data_header.IsUniform() == false &&
			   block.info_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
			{
				const int32_t global_id = block.footer.info_offsets[i].GetGlobalKey();
				//std::cerr << "info global: " << global_id << std::endl;
				if (this->memory_info[global_id].permanent_codec >= 0) {
					//std::cerr << "using permanent info" << std::endl;
					switch(this->memory_info[global_id].codec) {
						case(0): zstd_codec.Compress(block.info_containers[i]); break;
						case(1): this->CompressCodec1(block.info_containers[i]); break;
						case(2): this->CompressCodec2(block.info_containers[i]); break;
						case(3): this->CompressCodec3(block.info_containers[i]); break;
						case(4): this->CompressCodec4(block.info_containers[i]); break;
						case(5): this->CompressDelta(block.info_containers[i]); break;
					}
				} else if (this->memory_info[global_id].times_seen == 10) {
					//std::cerr << utility::timestamp("DEBUG") << "eval info " << global_id << std::endl;
					this->CompressEvaluate(block.info_containers[i], this->memory_info, global_id);
					this->memory_info[global_id].times_seen = 0;
				} else {
					//std::cerr << "reuse: " << (int)this->memory_info[i].codec << std::endl;
					switch(this->memory_info[global_id].codec) {
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

			if (i == YON_BLK_GT_INT16 || i == YON_BLK_GT_N_INT16) this->EncodeUnsignedVariableInt<uint16_t>(block.base_containers[i]);
			else if (i == YON_BLK_GT_INT32 || i == YON_BLK_GT_N_INT32) this->EncodeUnsignedVariableInt<uint32_t>(block.base_containers[i]);
			else
			*/
				zstd_codec.Compress(block.info_containers[i]);
		}
	}

	//io::BasicBuffer temp;
	for (uint32_t i = 0; i < block.footer.n_format_streams; ++i) {
		if (block.format_containers[i].header.n_entries) {
			if (block.format_containers[i].header.data_header.controller.type == YON_TYPE_FLOAT ||
			   block.format_containers[i].header.data_header.controller.type == YON_TYPE_DOUBLE) {
				zstd_codec.SetCompressionLevelData(general_level);
				zstd_codec.SetCompressionLevelStrides(general_level);
			}
			else {
				zstd_codec.SetCompressionLevel(general_level);
			}

			if (block.format_containers[i].header.n_entries) {
				/*
				const int32_t global_id = block.footer.format_offsets[i].GetGlobalKey();
				if (block.format_containers[i].header.data_header.IsUniform() == false &&
				   block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B)
				{
					if (this->memory_format[global_id].permanent_codec >= 0) {
						//std::cerr << "using permanent format" << std::endl;
						switch(this->memory_format[global_id].codec) {
							case(0): zstd_codec.Compress(block.format_containers[i]); break;
							case(1): this->CompressCodec1(block.format_containers[i]); break;
							case(2): this->CompressCodec2(block.format_containers[i]); break;
							case(3): this->CompressCodec3(block.format_containers[i]); break;
							case(4): this->CompressCodec4(block.format_containers[i]); break;
							case(5): this->CompressDelta(block.format_containers[i]); break;
						}
					} else if (this->memory_format[global_id].times_seen == 10) {
						//std::cerr << utility::timestamp("DEBUG") << "eval format " << global_id << std::endl;
						this->CompressEvaluate(block.format_containers[i], this->memory_format, global_id);
						this->memory_format[global_id].times_seen = 0;
					} else {
						//std::cerr << "reuse: " << (int)this->memory_format[i].codec << std::endl;
						switch(this->memory_format[global_id].codec) {
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
				if (block.format_containers[i].header.data_header.IsUniform() == false &&
				   (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B || block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_16B || block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_8B))
				{
					temp.resize(block.format_containers[i].data_uncompressed.size() + 65536);
					containers::StrideContainer s(block.format_containers[i]);

					uint8_t byte_size =

					for (int k = 0; k < s.size(); ++s) {
						block.format_containers[i].data_uncompressed
					}
				}
				*/


				/*
				if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_16B && block.format_containers[i].header.data_header.IsSigned() == false) {
					if (this->FormatStripedDelta<uint16_t, int32_t, uint32_t>(block.format_containers[i], n_samples) == false)
						this->EncodeUnsignedVariableInt<uint16_t>(block.format_containers[i]);


				} else if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B && block.format_containers[i].header.data_header.IsSigned() == false) {
					if (this->FormatStripedDelta<uint32_t, int64_t, uint64_t>(block.format_containers[i], n_samples) == false)
						this->EncodeUnsignedVariableInt<uint32_t>(block.format_containers[i]);

				}

				else if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_8B && block.format_containers[i].header.data_header.IsSigned() == false) {
					if (this->FormatStripedDelta<uint8_t, int16_t, uint16_t>(block.format_containers[i], n_samples) == false)
						zstd_codec.Compress(block.format_containers[i]);
				}


				else if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_32B && block.format_containers[i].header.data_header.IsSigned() == true) {
					this->EncodeZigZagVariableInt32(block.format_containers[i]);

				} else if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_16B && block.format_containers[i].header.data_header.IsSigned() == true) {
					if (this->FormatStripedDelta<int16_t, int32_t, uint32_t>(block.format_containers[i], n_samples) == false)
						this->EncodeZigZagVariableInt16(block.format_containers[i]);
					//	this->EncodeZigZag<int16_t, uint16_t>(block.format_containers[i]);

				} else if (block.format_containers[i].header.data_header.GetPrimitiveType() == YON_TYPE_8B && block.format_containers[i].header.data_header.IsSigned() == true) {
					//
					if (this->FormatStripedDelta<int8_t, int16_t, uint16_t>(block.format_containers[i], n_samples) == false)
						this->EncodeZigZag<int8_t, uint8_t>(block.format_containers[i]);

				}
				else
				*/
					zstd_codec.Compress(block.format_containers[i]);

			}
		}
	}
	return true;
}

bool CompressionManager::Decompress(variant_block_type& block) {
	if (block.base_containers[YON_BLK_PPA].GetSizeCompressed() &&
	   block.base_containers[YON_BLK_PPA].header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
	{
		if (!this->Decompress(block.base_containers[YON_BLK_PPA], *block.gt_ppa)) {
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress GT permutation information!" << std::endl;
			return false;
		}
	}

	for (uint32_t i = 1; i < YON_BLK_N_STATIC; ++i) {
		if (block.base_containers[i].GetSizeCompressed() &&
		   block.base_containers[i].IsEncrypted() == false)
		{
			if (!this->Decompress(block.base_containers[i])) {
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress basic container!" << std::endl;
				return false;
			}
		}
	}

	for (uint32_t i = 0; i < block.footer.n_info_streams; ++i) {
		if (block.info_containers[i].GetSizeCompressed() &&
		   block.info_containers[i].IsEncrypted() == false)
		{
			if (!this->Decompress(block.info_containers[i])) {
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress INFO container " << i << "/" << block.footer.n_info_streams << "!" << std::endl;
				return false;
			}
		}
	}

	for (uint32_t i = 0; i < block.footer.n_format_streams; ++i) {
		if (block.format_containers[i].GetSizeCompressed() &&
		   block.format_containers[i].IsEncrypted() == false)
		{
			if (!this->Decompress(block.format_containers[i])) {
				std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress FORMAT container " << i << "/" << block.footer.n_format_streams << "!" << std::endl;
				return false;
			}
		}
	}

	return true;
}

bool CompressionManager::Decompress(container_type& container, yon_gt_ppa& gt_ppa) {
	return(this->zstd_codec.Decompress(container, gt_ppa));
}

bool CompressionManager::Decompress(container_type& container) {
	// Ascertain that data is not encrypted
	if (container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE) {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if (container.header.data_header.controller.encoder == YON_ENCODE_ZSTD) {
		if (!this->zstd_codec.Decompress(container)) {
			std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress data!" << std::endl;
			return false;
		}
	} else if (container.header.data_header.controller.encoder == YON_ENCODE_NONE) {
		if (!this->no_codec.Decompress(container)) {
			std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress data!" << std::endl;
			return false;
		}
	} else if (container.header.data_header.controller.encoder == YON_ENCODE_ZPAQ) {
		std::cerr << utility::timestamp("ERROR","CODEC-ZPAQ") << "ZPAQ is no longer supported!" << std::endl;
		return false;
	} else {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to decompress! Illegal codec!" << std::endl;
		return false;
	}

	if (container.header.data_header.controller.mixedStride) {
		if (container.header.stride_header.controller.encoder == YON_ENCODE_ZSTD) {
			if (!this->zstd_codec.DecompressStrides(container)) { std::cerr << utility::timestamp("ERROR","CODEC-ZSTD") << "Failed to decompress strides!" << std::endl; return false; }
		} else if (container.header.stride_header.controller.encoder == YON_ENCODE_NONE) {
			if (!this->no_codec.DecompressStrides(container)) { std::cerr << utility::timestamp("ERROR","CODEC-NONE") << "Failed to decompress strides!" << std::endl; return false; }
		} else if (container.header.stride_header.controller.encoder == YON_ENCODE_ZPAQ) {
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
