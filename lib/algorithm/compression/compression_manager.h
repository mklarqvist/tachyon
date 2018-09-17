#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "containers/variant_block.h"
#include "containers/stride_container.h"
#include "algorithm/compression/fastdelta.h"

namespace tachyon{
namespace algorithm{

class CompressionManager{
public:
	typedef CompressionManager        self_type;
	typedef UncompressedCodec         no_codec_type;
	typedef ZSTDCodec                 zstd_codec_type;
	typedef containers::VariantBlock  variant_block_type;
	typedef containers::DataContainer container_type;
	typedef containers::VariantBlockFooter footer_type;

public:
	CompressionManager() = default;
	~CompressionManager() = default;

	bool Compress(variant_block_type& block, const uint8_t general_level, const uint32_t n_samples);
	bool Decompress(variant_block_type& block);
	bool Decompress(container_type& container, yon_gt_ppa& gt_ppa);

	template <class int_t>
	bool FormatStrideShift(container_type& container, const uint32_t n_samples){
		containers::StrideContainer<> s(container);
		const int_t* src = reinterpret_cast<const int_t*>(container.data_uncompressed.data());
		const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int_t);
		assert(container.data_uncompressed.size() % sizeof(int_t) == 0);

		io::BasicBuffer dst_buffer(n_entries*sizeof(int_t) + 65536);
		int_t* dst = reinterpret_cast<int_t*>(dst_buffer.data());

		uint32_t offset = 0;

		int_t** buckets = new int_t*[100];

		// Iterate over available entries.
		for(int i = 0; i < s.size(); ++i){
			const uint32_t n_width  = n_samples*s[i]; // n_entries in row
			const int_t* l_src = &src[offset]; // current src data
			int_t* l_dst = &dst[offset]; // current src data

			// Setup buckets
			for(int b = 0; b < s[i]; ++b){
				buckets[b] = &l_dst[n_samples*b];
			}

			/*
			if(s[i] == 3){
				std::cerr << "stride: " << s[i] << std::endl;
				for(int k = 0; k < n_width; ++k){
					std::cerr << "," << (int)l_src[k];
				}
				std::cerr << std::endl;
			}
			*/

			// First
			uint32_t l_offset = 0;

			// Iterate over stride
			for(int j = 0; j < n_samples; ++j){
				// Iterate over current stride.
				for(int b = 0; b < s[i]; ++b){
					*buckets[b] = l_src[l_offset++];
					++buckets[b];
				}
			}

			/*
			if(s[i] == 3){
				for(int k = 0; k < n_width; ++k){
					std::cerr << "," << (int)l_dst[k];
				}
				std::cerr << std::endl;
				exit(1);
			}
			*/

			assert(l_offset == n_width);
			offset += n_width;
		}

		assert(offset == n_entries);
		const uint32_t b_original = container.data_uncompressed.size();
		dst_buffer.n_chars_ = container.data_uncompressed.size();

		delete [] buckets;

		/*
		std::cout << "original:" << std::endl;
		for(int i = 0; i < n_entries; ++i){
			std::cout << "," << (int)src[i];
		}
		std::cout << std::endl;
		std::cout << "delta:" << std::endl;
		for(int i = 0; i < n_entries; ++i){
			std::cout << "," << (int)dst[i];
		}
		std::cout << std::endl;
		exit(1);
		*/

		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-pack-" << sizeof(int_t)*8 << "-default: " <<  container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;
		const uint32_t std_size = container.data.size();

		io::BasicBuffer ubackup; ubackup = std::move(container.data_uncompressed);
		//io::BasicBuffer cbackup; cbackup = std::move(container.data);

		container.data_uncompressed = std::move(dst_buffer);
		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-pack-" << sizeof(int_t)*8 << ": " << b_original << "->" << container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;
		return true;
	}

	template <class in_int_t, class out_int_t, class out_uint_t>
	bool FormatHorizontalDelta(container_type& container, const uint32_t n_samples){
		containers::StrideContainer<> s(container);
		const in_int_t* src = reinterpret_cast<const in_int_t*>(container.data_uncompressed.data());
		const uint32_t n_entries = container.data_uncompressed.size() / sizeof(in_int_t);
		assert(container.data_uncompressed.size() % sizeof(in_int_t) == 0);

		io::BasicBuffer dst_buffer(n_entries*sizeof(out_uint_t) + 65536);
		out_uint_t* dst = reinterpret_cast<out_uint_t*>(dst_buffer.data());

		uint32_t offset = 0;
		const uint8_t shift = sizeof(out_int_t)*8 - 1;

		// Iterate over available entries.
		uint32_t last_stride = s[0];
		uint32_t skips = 0;
		const uint32_t ref_stride = s[0];
		for(int i = 0; i < s.size(); ++i){
			bool perform_check = false;
			if(s[i] == ref_stride){
				perform_check = true;
			} else {
				skips += n_samples*s[i];
			}

			const uint32_t n_width  = n_samples*s[i]; // n_entries in row
			const in_int_t* l_src = &src[offset]; // current src data
			out_uint_t* l_dst = &dst[offset]; // current src data

			if(perform_check && i != 0){
				//std::cerr << "perfomring check" << std::endl;
				// Iterate over stride
				uint32_t l_offset = 0;
				for(int j = 0; j < n_samples; ++j){
					// Iterate over current stride.
					for(int bucket = 0; bucket < s[i]; ++bucket){
						// Compute prefix sum by looking at n_samples*s[i] into the past.
						//std::cerr << "check: " << (int)src[offset + l_offset] << " and " << (int)src[((int)offset + l_offset) - n_width] << "=" << (int)(out_int_t)src[offset + l_offset] - (out_int_t)src[(offset + l_offset) - n_width] << std::endl;
						out_int_t delta = (out_int_t)src[offset + l_offset] - (out_int_t)src[(offset + l_offset) - (n_width + skips)];
						l_dst[l_offset] = out_uint_t((delta << 1) ^ (delta >> shift));
						++l_offset;
					}
				}
				assert(l_offset == n_width);
				skips = 0;

			} else {
				//std::cerr << "no checks done copying data: " << s[i] << "/" << last_stride << "@" << i << std::endl;
				for(int j = 0; j < n_width; ++j)
					l_dst[j] = l_src[j];
			}
			last_stride = s[i];
			offset += n_width;
		}

		assert(offset == n_entries);
		const uint32_t b_original = container.data_uncompressed.size();
		dst_buffer.n_chars_ = n_entries*sizeof(out_int_t);
		assert(dst_buffer.n_chars_ > container.data_uncompressed.size());

		/*
		std::cout << "original:" << std::endl;
		for(int i = 0; i < n_entries; ++i){
			std::cout << "," << (int)src[i];
		}
		std::cout << std::endl;
		std::cout << "delta:" << std::endl;
		for(int i = 0; i < n_entries; ++i){
			std::cout << "," << (int)dst[i];
		}
		std::cout << std::endl;
		exit(1);
		*/

		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-fDdelta-" << sizeof(in_int_t)*8 << "-default: " <<  container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;
		const uint32_t std_size = container.data.size();

		io::BasicBuffer ubackup; ubackup = std::move(container.data_uncompressed);
		//io::BasicBuffer cbackup; cbackup = std::move(container.data);

		container.data_uncompressed = std::move(dst_buffer);
		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-fDdelta-" << sizeof(in_int_t)*8 << ": " << b_original << "->" << container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;
		return true;
	}

	template <class in_int_t, class out_int_t, class out_uint_t>
	bool FormatStripedDelta(container_type& container, const uint32_t n_samples){
		containers::StrideContainer<> s(container);
		const in_int_t* src = reinterpret_cast<const in_int_t*>(container.data_uncompressed.data());
		const uint32_t n_entries = container.data_uncompressed.size() / sizeof(in_int_t);
		assert(container.data_uncompressed.size() % sizeof(in_int_t) == 0);

		io::BasicBuffer dst_buffer(n_entries*sizeof(out_uint_t) + 65536);
		out_uint_t* dst = reinterpret_cast<out_uint_t*>(dst_buffer.data());

		uint32_t offset = 0;
		const uint8_t shift = sizeof(out_int_t)*8 - 1;
		// Iterate over available entries.
		for(int i = 0; i < s.size(); ++i){
			const uint32_t n_width  = n_samples*s[i]; // n_entries in row
			const in_int_t* l_src = &src[offset]; // current src data
			out_uint_t* l_dst = &dst[offset]; // current src data

			// First
			uint32_t l_offset = 0;
			for(int stride = 0; stride < s[i]; ++stride){
				l_dst[l_offset] = l_src[l_offset];
				++l_offset;
			}

			// Iterate over stride
			for(int j = 1; j < n_samples; ++j){
				// Iterate over current stride.
				for(int bucket = 0; bucket < s[i]; ++bucket){
					// Compute prefix sum.
					out_int_t delta = (out_int_t)l_src[l_offset] - (out_int_t)l_src[l_offset - s[i]];
					l_dst[l_offset] = out_uint_t((delta << 1) ^ (delta >> shift));
					++l_offset;
				}
			}

			assert(l_offset == n_width);
			offset += n_width;
		}
		assert(offset == n_entries);
		const uint32_t b_original = container.data_uncompressed.size();
		dst_buffer.n_chars_ = n_entries*sizeof(out_int_t);
		//assert(dst_buffer.n_chars_ > container.data_uncompressed.size());

		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-fdelta-" << sizeof(in_int_t)*8 << "-default: " <<  container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;
		const uint32_t std_size = container.data.size();

		io::BasicBuffer ubackup; ubackup = std::move(container.data_uncompressed);
		//io::BasicBuffer cbackup; cbackup = std::move(container.data);

		container.data_uncompressed = std::move(dst_buffer);
		zstd_codec.Compress(container);
		std::cerr << utility::timestamp("DEBUG") << "zstd-fdelta-" << sizeof(in_int_t)*8 << ": " << b_original << "->" << container.data_uncompressed.size() << "->" << container.data.size() << " -> " << (float)b_original/container.data.size() << std::endl;

		if((float)b_original/container.data.size() < 5){
			if(sizeof(in_int_t) == 1){
				if(container.data.size() < std_size)
					return true;
			}
			std::cerr << utility::timestamp("DEBUG") << "Attempt backup..." << std::endl;
			container.data_uncompressed = std::move(ubackup);
			return false;
		} else {
			std::cerr << utility::timestamp("LOG","COMPRESSION") << "Delta: " << (float)b_original/container.data.size() << std::endl;
		}


		return true;
	}

	bool EncodeZigZagVariableInt32(container_type& container){
		if(container.header.data_header.GetPrimitiveType() == YON_TYPE_32B){
			io::BasicBuffer buf(container.data_uncompressed.size() + 65536);
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

			io::BasicBuffer restore(container.data_uncompressed.size() + 65536);
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

			//io::BasicBuffer backup;
			//backup = std::move(container.data_uncompressed);

			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG") << "zstd-vintz32: " << b_original << "->" << container.data.size() << "->" << (float)b_original/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);
		}
		return true;
	}

	bool EncodeZigZagVariableInt16(container_type& container){
		if(container.header.data_header.GetPrimitiveType() == YON_TYPE_16B){
			io::BasicBuffer buf(container.data_uncompressed.size() + 65536);
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

			io::BasicBuffer restore(container.data_uncompressed.size() + 65536);
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

			//io::BasicBuffer backup;
			//backup = std::move(container.data_uncompressed);
			container.data_uncompressed = std::move(buf);
			container.header.data_header.uLength = buf.size();
			container.header.data_header.controller.preprocessor = 1;
			zstd_codec.Compress(container);
			std::cerr << utility::timestamp("DEBUG")<< "zstd-vintz16 " << b_original << "->" << container.data.size() << "->" << (float)b_original/container.data.size() << std::endl;
			//container.data_uncompressed = std::move(backup);
		}
		return true;
	}

	template <class int_t, class uint_t>
	bool EncodeZigZag(container_type& container){
			io::BasicBuffer buf(container.data_uncompressed.size() + 65536);
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

			//io::BasicBuffer backup;
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
			io::BasicBuffer buf(container.data_uncompressed.size() + 65536);
			const uint16_t* data = reinterpret_cast<const uint16_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int16_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			for(int i = 0; i < n_entries; ++i)
				buf += (data[i] << 1) ^ (data[i] >> 15);

			buf.n_chars_ = container.data_uncompressed.size();


			//zstd_codec.Compress(container);
			//std::cerr << "zstd-ref-zigzag8o: " << container.data.size() << " " << (float)container.data_uncompressed.size()/container.data.size() << std::endl;

			//io::BasicBuffer backup;
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

		io::BasicBuffer buf(container.data_uncompressed.size() + 65536);
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

		//io::BasicBuffer backup;
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
