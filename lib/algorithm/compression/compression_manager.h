#ifndef ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSION_MANAGER_H_

#include "uncompressed_codec.h"
#include "zstd_codec.h"
#include "variant_container.h"
#include "algorithm/compression/compression_manager.h"

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

	bool EncodeZigZagVariableInt32(container_type& container) {
		if (container.header.data_header.GetPrimitiveType() == YON_TYPE_32B) {
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint32_t* data = reinterpret_cast<const uint32_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int32_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			uint32_t bytes = 0;
			for (int i = 0; i < n_entries; ++i) {
				// No benefit. Fallback.
				if (bytes + 10 >= buf.capacity()) {
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
			for (int i = 0; i < n_entries; ++i) {
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
		}
		return true;
	}

	bool EncodeZigZagVariableInt16(container_type& container) {
		if (container.header.data_header.GetPrimitiveType() == YON_TYPE_16B) {
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint16_t* data = reinterpret_cast<const uint16_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int16_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			uint32_t bytes = 0;
			for (int i = 0; i < n_entries; ++i) {
				// No benefit. Fallback.
				if (bytes + 10 >= buf.capacity()) {
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
			for (int i = 0; i < n_entries; ++i) {
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
		}
		return true;
	}

	template <class int_t, class uint_t>
	bool EncodeZigZag(container_type& container) {
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const int_t* data = reinterpret_cast<const int_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			const uint8_t shift = sizeof(int_t)*8 - 1;
			for (int i = 0; i < n_entries; ++i) {
				// If debugging:
				/*
				uint_t a = ((data[i] << 1) ^ (data[i] >> shift));
				int_t  b = (a >> 1) ^ -(a & 1);
				if (data[i] != b) {
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

	bool EncodeZigZagInt16(container_type& container) {
		if (container.header.data_header.GetPrimitiveType() == YON_TYPE_16B) {
			yon_buffer_t buf(container.data_uncompressed.size() + 65536);
			const uint16_t* data = reinterpret_cast<const uint16_t*>(container.data_uncompressed.data());
			const uint32_t n_entries = container.data_uncompressed.size() / sizeof(int16_t);
			uint8_t* out = reinterpret_cast<uint8_t*>(buf.data());

			for (int i = 0; i < n_entries; ++i)
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
	bool EncodeUnsignedVariableInt(container_type& container) {
		if (container.data_uncompressed.size() == 0) {
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
		for (int i = 0; i < n_entries; ++i) {
			// No benefit. Fallback.
			if (bytes + 10 >= buf.capacity()) {
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
