#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "variant_container.h"
#include "algorithm/compression/compression_manager.h"

#include "algorithm/compression/packed_array.h"

namespace tachyon{
namespace algorithm{

class CompressionManager{
public:
	typedef CompressionManager self_type;
	typedef UncompressedCodec  no_codec_type;
	typedef ZSTDCodec          zstd_codec_type;
	typedef yon1_vb_t          variant_block_type;
	typedef yon1_dc_t          container_type;
	typedef yon_vb_ftr         footer_type;

public:
	CompressionManager() = default;
	~CompressionManager() = default;

	bool Compress(variant_block_type& block, const uint8_t general_level, const uint32_t n_samples);
	bool Decompress(variant_block_type& block);
	bool Decompress(container_type& container, yon_gt_ppa& gt_ppa);

	bool EncodeHorizontalDelta(container_type& container, const uint32_t n_samples){
		if(container.header.data_header.IsSigned() == false){
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B): return(this->EncodeHorizontalDelta_<uint8_t>(container,n_samples));
			case(YON_TYPE_16B): return(this->EncodeHorizontalDelta_<uint16_t>(container,n_samples));
			case(YON_TYPE_32B): return(this->EncodeHorizontalDelta_<uint32_t>(container,n_samples));
			case(YON_TYPE_64B): return(this->EncodeHorizontalDelta_<uint64_t>(container,n_samples));
			default: return false;
			}
		}
		else {
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B): return(this->EncodeHorizontalDelta_<int8_t>(container,n_samples));
			case(YON_TYPE_16B): return(this->EncodeHorizontalDelta_<int16_t>(container,n_samples));
			case(YON_TYPE_32B): return(this->EncodeHorizontalDelta_<int32_t>(container,n_samples));
			case(YON_TYPE_64B): return(this->EncodeHorizontalDelta_<int64_t>(container,n_samples));
			default: return false;
			}
		}

		return false;
	}

	template <class int_t>
	bool EncodeHorizontalDelta_(container_type& container, const uint32_t n_samples){
		if(container.data_uncompressed.size() == 0) return false;

		const uint32_t nrdata = container.data_uncompressed.size() / sizeof(int_t);
		assert(container.data_uncompressed.size() % sizeof(int_t) == 0);

		if(container.header.data_header.HasMixedStride() == false){
			std::cerr << "no mxied stride" << std::endl;
			return false;
			//exit(1);
		} else {
			yon_cont_ref_iface* it = nullptr;
			switch(container.header.stride_header.controller.type){
			case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
			case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
			case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
			case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
			}
			assert(it != nullptr);

			uint32_t n_stride[256]; memset(n_stride, 0, sizeof(uint32_t)*256);
			uint32_t n_largest_stride = 0;
			for(int i = 0; i < it->n_elements_; ++i){
				//std::cerr << "stride-" << it->GetInt32(i) << std::endl;
				const uint32_t cur_stride = it->GetInt32(i);
				n_largest_stride = cur_stride > n_largest_stride ? cur_stride : n_largest_stride;
				for(int j = 0; j < it->GetInt32(i); ++j)
					++n_stride[j];
			}
			++n_largest_stride;

			std::cerr << "input=" << container.data_uncompressed.size() << std::endl;

			const int_t* rdata  = reinterpret_cast<const int_t*>(container.data_uncompressed.data());
			yon1_dc_t* streams  = new yon1_dc_t[n_largest_stride];
			uint64_t* s_offsets = new uint64_t[n_largest_stride];
			int_t** tdata       = new int_t*[n_largest_stride];
			for(int i = 0; i < n_largest_stride; ++i){
				std::cerr << "stride=" << n_stride[i] << std::endl;
				streams[i].data_uncompressed.resize(n_samples * n_stride[i] * sizeof(int_t));
				std::cerr << "capacity=" << streams[i].data_uncompressed.capacity() << std::endl;
				tdata[i] = reinterpret_cast<int_t*>(streams[i].data_uncompressed.data());
			}
			memset(s_offsets, 0, sizeof(uint64_t)*n_largest_stride);

			uint32_t offset = 0;
			for(int i = 0; i < it->n_elements_; ++i){
				uint32_t local_offset = 0;
				const int_t* ldata = &rdata[offset];
				uint32_t cur_stride = it->GetInt32(i);

				for(int s = 0; s < n_samples; ++s){ // over samples
					const uint32_t cur_stride = it->GetInt32(i);
					for(int j = 0; j < cur_stride; ++j, ++local_offset){ // over stride
						tdata[j][s_offsets[j]++] = ldata[local_offset];
						++tdata[j][s_offsets[j]++];
						//if(j == 1) std::cerr << "adding to=" << j << "/" << cur_stride << " at " << s_offsets[j] << " value=" << (int64_t)ldata[local_offset] << std::endl;
					}
				}
				offset += sizeof(int_t)*n_samples*cur_stride;
			}
			assert(offset == container.data_uncompressed.size());

			uint32_t n_total = 0;
			for(int i = 0; i < n_largest_stride; ++i){
				if(s_offsets[i] == 0) continue;
				streams[i].data_uncompressed.n_chars_ = s_offsets[i]*sizeof(int_t);
				std::cerr << "checking=" << streams[i].data_uncompressed.size() << "/" << streams[i].data_uncompressed.capacity() << std::endl;

				//streams[i].UpdateContainer(true,true);
				this->zstd_codec.Compress(streams[i]);
				std::cerr << streams[i].data_uncompressed.size() << "->" << streams[i].data.size() << "\t" << (float)streams[i].data_uncompressed.size()/streams[i].data.size() << std::endl;

				if((float)streams[i].data_uncompressed.size()/streams[i].data.size() < 5){
					int_t maxv = std::numeric_limits<int_t>::min(), minv = std::numeric_limits<int_t>::max(),
							maxv2 = std::numeric_limits<int_t>::min(), minv2 = std::numeric_limits<int_t>::max();

					/*
					uint32_t nbins = s_offsets[i] / 512;
					uint32_t o_tot = 0;
					for(int j = 0; j < nbins; ++j){
						for(int k = 0; k < 512; ++k, ++o_tot){
							maxv = tdata[i][o_tot] > maxv ? tdata[i][o_tot] : maxv;
							minv = tdata[i][o_tot] < minv ? tdata[i][o_tot] : minv;
							minv2 = tdata[i][o_tot] < minv2 && tdata[i][o_tot] != minv ? tdata[i][o_tot] : minv2;
							maxv2 = tdata[i][o_tot] > maxv2 && tdata[i][o_tot] != maxv ? tdata[i][o_tot] : maxv2;
						}
						//std::cerr << "range=" << (int64_t)minv << "-" << (int64_t)maxv << " 2nd-min=" << (int64_t)minv2 << " 2nd-max=" << (int64_t)maxv2 << " " << log2((maxv-minv)+1) << std::endl;
						maxv = std::numeric_limits<int_t>::min(), minv = std::numeric_limits<int_t>::max();
						maxv2 = std::numeric_limits<int_t>::min(), minv2 = std::numeric_limits<int_t>::max();
					}
					*/

					//if(i == 1){
					for(int j = 0; j < s_offsets[i]; ++j){
						maxv = tdata[i][j] > maxv ? tdata[i][j] : maxv;
						minv = tdata[i][j] < minv && tdata[i][j] > std::numeric_limits<int_t>::min()+1 ? tdata[i][j] : minv;
						//std::cerr << (int)tdata[i][j] << ",";
					}
					//std::cerr << std::endl;
					std::cerr << "range=" << (int64_t)minv << "-" << (int64_t)maxv << " at " << sizeof(int_t) << " sign=" << container.header.data_header.IsSigned() << std::endl;

					if(container.header.data_header.IsSigned() == false){
						if(sizeof(int_t) > 1)
							this->EncodeUnsignedVariableInt<int_t>(streams[i]);
					}

					//}
					//exit(1);
				}

				n_total += streams[i].data.size();
				//delete tdata[i];
			}
			std::cerr << "total: " << container.data_uncompressed.size() << "->" << n_total << "\t" << (float)container.data_uncompressed.size()/n_total << std::endl;
			//std::cerr << "bitpack: " << a->count << std::endl;


			this->zstd_codec.Compress(container);
			std::cerr << "ref: " << container.data_uncompressed.size() << "->" << container.data.size() << "\t" << (float)container.data_uncompressed.size()/container.data.size() << std::endl;

			//delete a;
			delete [] tdata;
			delete [] s_offsets;
			delete [] streams;
			delete it;
			return true;
		}
	}

	bool EncodeZigZagVariableInt32(container_type& container){
		//if(container.header.data_header.GetPrimitiveType() == YON_TYPE_32B){
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint32_t* data = reinterpret_cast<const uint32_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int32_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			uint32_t bytes = 0;
			for(int i = 0; i < n_entries; ++i){
				// No benefit. Fallback.
				if(bytes + 10 >= buf.capacity()){
					std::cerr << utility::timestamp("DEBUG") << "Fallback in vint32z: " << bytes << ">" << buf.capacity() << std::endl;
					zstd_codec.Compress(container);
					return true;
				}
				bytes += EncodeVarint<uint32_t>((data[i] << 1) ^ (data[i] >> 31), &out[bytes]);
			}

			buf.n_chars_ = bytes;
			const uint32_t b_original = container.data_uncompressed.size();
			/*
			std::cerr << "bytes: " << container.data_uncompressed.size() << "->" << bytes << std::endl;

			const uint64_t hash_orignal = XXH64(container.data_uncompressed.data(), container.data_uncompressed.size(), 0);

			yon_buffer_t restore(container.data_uncompressed.size() + 65536);
			size_t offset = 0;
			for(int i = 0; i < n_entries; ++i){
				restore += DecodeZigzag16(DecodeVarint<uint16_t>(out, offset));
			}
			const int16_t* data_rst = reinterpret_cast<const int16_t*>(restore.data());
			//std::cerr << "restore: " << restore.size() << "==" << container.data_uncompressed.size() << std::endl;
			const uint64_t hash_restore = XXH64(restore.data(), restore.size(), 0);
			std::cerr << "hash comp: " << hash_orignal << "==" << hash_restore << std::endl;

			assert(hash_orignal == hash_restore);
*/

			//yon_buffer_t backup;
			//backup = std::move(container.data_uncompressed);

			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG") << "zstd-vintz32: " << b_original << "->" << container.data.size() << "->" << (float)b_original/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);
		//}
		return true;
	}

	bool EncodeZigZagVariableInt16(container_type& container){
		//if(container.header.data_header.GetPrimitiveType() == YON_TYPE_16B){
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint16_t* data = reinterpret_cast<const uint16_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int16_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			uint32_t bytes = 0;
			for(int i = 0; i < n_entries; ++i){
				// No benefit. Fallback.
				if(bytes + 10 >= buf.capacity()){
					std::cerr << utility::timestamp("DEBUG") << "Fallback in vint16z: " << bytes << ">" << buf.capacity() << std::endl;
					zstd_codec.Compress(container);
					return true;
				}
				bytes += EncodeVarint<uint16_t>((data[i] << 1) ^ (data[i] >> 15), &out[bytes]);
			}
			buf.n_chars_ = bytes;
			const uint32_t b_original = container.data_uncompressed.size();
			/*
			std::cerr << "bytes: " << container.data_uncompressed.size() << "->" << bytes << std::endl;

			const uint64_t hash_orignal = XXH64(container.data_uncompressed.data(), container.data_uncompressed.size(), 0);

			yon_buffer_t restore(container.data_uncompressed.size() + 65536);
			size_t offset = 0;
			for(int i = 0; i < n_entries; ++i){
				restore += DecodeZigzag16(DecodeVarint<uint16_t>(out, offset));
			}
			const int16_t* data_rst = reinterpret_cast<const int16_t*>(restore.data());
			//std::cerr << "restore: " << restore.size() << "==" << container.data_uncompressed.size() << std::endl;
			const uint64_t hash_restore = XXH64(restore.data(), restore.size(), 0);
			std::cerr << "hash comp: " << hash_orignal << "==" << hash_restore << std::endl;

			assert(hash_orignal == hash_restore);
*/

			//yon_buffer_t backup;
			//backup = std::move(container.data_uncompressed);
			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG")<< "zstd-vintz16 " << b_original << "->" << container.data.size() << "->" << (float)b_original/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);
		//}
		return true;
	}

	template <class int_t, class uint_t>
	bool EncodeZigZag(container_type& container){
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const int_t* data = reinterpret_cast<const int_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			const uint8_t shift = sizeof(int_t)*8 - 1;
			for(int i = 0; i < n_entries; ++i){
				// If debugging:
				/*
				uint_t a = ((data[i] << 1) ^ (data[i] >> shift));
				int_t  b = (a >> 1) ^ -(a & 1);
				if(data[i] != b){
					std::cerr << i << " Error: input=" << (int)data[i] << "!=" << (int)b << std::endl;
					exit(1);
				}
				*/
				buf += (uint_t)((data[i] << 1) ^ (data[i] >> shift));
			}

			buf.n_chars_ = container.data_uncompressed.size();


			//zstd_codec.Compress(container);
			//std::cerr << "zstd-ref-zigzag8o: " << container.data.size() << " " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;

			//yon_buffer_t backup;
			//backup = std::move(container.data_uncompressed);
			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG") << "zstd-zig" << sizeof(int_t)*8 << "o: " << container.data_uncompressed.size() << "->" << container.data.size() << " " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);

		return true;
	}

	bool EncodeZigZagInt16(container_type& container){
		if(container.header.data_header.GetPrimitiveType() == YON_TYPE_16B){
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint16_t* data = reinterpret_cast<const uint16_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int16_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			for(int i = 0; i < n_entries; ++i)
				buf += (data[i] << 1) ^ (data[i] >> 15);

			buf.n_chars_ = container.data_uncompressed.size();


			//zstd_codec.Compress(container);
			//std::cerr << "zstd-ref-zigzag8o: " << container.data.size() << " " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;

			//yon_buffer_t backup;
			//backup = std::move(container.data_uncompressed);
			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG") << "zstd-zig16o: " << container.data_uncompressed.size() << "->" << container.data.size() << " " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);
		}
		return true;
	}

	template <class int_t>
	bool EncodeUnsignedVariableInt(container_type& container){
		if(container.data_uncompressed.size() == 0){
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.data.n_chars_                         = 0;
			container.header.data_header.cLength            = 0;
			container.header.data_header.uLength            = 0;
			return false;
		}

		//std::cerr << utility::timestamp("DEBUG") << "VarInt" << sizeof(int_t)*8 << std::endl;
		zstd_codec.Compress(container);
		const uint32_t base = container.data.size();

		yon_buffer_t buf(container.data_uncompressed.size() + 65536);
		const int_t* data = reinterpret_cast<const int_t*>(container.data_uncompressed.data());
		const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int_t);
		uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

		uint32_t bytes = 0;
		for(int i = 0; i < n_entries; ++i){
			// No benefit. Fallback.
			if(bytes + 10 >= buf.capacity()){
				zstd_codec.Compress(container);
				std::cerr << utility::timestamp("DEBUG") << "Fallback in uvint" << sizeof(int_t)*8 << ": " << bytes << ">" << buf.capacity() << "-> fallback: " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;

				return true;
			}
			bytes += EncodeVarint<int_t>(data[i], &out[bytes]);
		}
		//std::cerr << "bytes: " << container.data_uncompressed.size() << "->" << bytes << std::endl;
		buf.n_chars_ = bytes;
		const uint32_t b_original = container.data_uncompressed.size();

		//yon_buffer_t backup;
		//backup = std::move(container.data_uncompressed);
		container.data_uncompressed = std::move(buf);
		container.header.data_header.uLength = buf.size();
		container.header.data_header.controller.preprocessor = 1;
		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "varint-" << sizeof(int_t)*8 << ": " << b_original << "->" << container.data.size() << "->" << (float)b_original/container.data.size() << " default is: " << (float)b_original/base << std::endl;
		//container.data_uncompressed = std::move(backup);
		return true;
	}

	/**<
	 * Decompress a data container.
	 * @param container Src and dst target container.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool Decompress(container_type& container);

public:
	no_codec_type   no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_ */
