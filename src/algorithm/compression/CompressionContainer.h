#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"
#include "zstd.h"
#include "common/zstd_errors.h"
#include "../../third_party/r16N.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"

namespace Tachyon{
namespace Compression{

#define MIN_COMPRESSION_FOLD 1.05

inline bool bytePreprocessBits(const char* const data, const U32& size, char* destination){
	if(size == 0) return false;

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data);
	BYTE* target[8];
	const U32 s = size/8;
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &dest[s*i];

	U32 k = 0; U32 p = 0;
	for(U32 i = 0; i < size; ++i, ++k){
		for(U32 j = 0; j < 8; ++j)
			target[j][p] |= ((d[i] & (1 << j)) >> j) << k;

		if(i % 7 == 0){ k = 0; ++p; }
	}

	return true;
}

class CompressionContainer{
private:
	typedef CompressionContainer self_type;
protected:
	typedef Core::StreamContainer stream_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::PermutationManager permutation_type;

public:
	CompressionContainer(){}
	virtual ~CompressionContainer(){}
	virtual const bool encode(permutation_type& manager) =0;
	virtual const bool encode(stream_type& stream) =0;
	virtual const bool encodeStrides(stream_type& stream) =0;
	virtual const bool decode(stream_type& stream) =0;

protected:
	buffer_type buffer;
};

// GZIP of delfate codec
class ZSTDCodec : public CompressionContainer{
private:
	typedef ZSTDCodec self_type;
	typedef IO::TGZFController controller_type;

public:
	ZSTDCodec() : compression_level(0){}
	~ZSTDCodec(){ }

	void setCompressionLevel(const int& c){ this->compression_level = c; }

	const bool assess(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer + 65536);
		//for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_strides.data, stream.buffer_strides.pointer, this->compression_level);
			if(ZSTD_isError(ret)){
				std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
				//exit(1);
			}
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << this->compression_level << ": " << stream.buffer_strides.pointer << '\t' << ret << '\t' << (double)stream.buffer_strides.pointer/ret << "-fold" << std::endl;
		//}

		U32 out_size = this->buffer.capacity();
		rans_compress_to_4x16((BYTE*)stream.buffer_strides.data, stream.buffer_strides.pointer, (BYTE*)this->buffer.data, &out_size, 0);
		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "RANS0\t" << stream.buffer_strides.pointer << "\t" << out_size << '\t' << (double)stream.buffer_strides.pointer/out_size << std::endl;

		/*
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)stream.buffer_data.data, stream.buffer_data.pointer);



		buffer_type test(this->buffer);
		U32 out2 = stream.buffer_data.pointer;
		rans_uncompress_to_4x16((BYTE*)this->buffer.data, out_size, (BYTE*)&test.data[0], &out2, 0);
		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "RANS0-RE\t" << out_size << "\t" << out2 << std::endl;

		U32 crc2 = crc32(0, NULL, 0);
		crc2 = crc32(crc2, (Bytef*)test.data, out2);

		std::cerr << crc << '\t' << crc2 << '\t' << (crc == crc2) << std::endl;

		test.deleteAll();
	*/

		/*
		out_size = this->buffer.capacity();
		rans_compress_to_4x16((BYTE*)&stream.buffer_data.data[stream.buffer_data.pointer-1], stream.buffer_data.pointer, (BYTE*)this->buffer.data, &out_size, 1);
		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "RANS1\t" << stream.buffer_data.pointer << "\t" << out_size << '\t' << (double)stream.buffer_data.pointer/out_size << std::endl;
		*/
		return true;
	}

	const bool encode(stream_type& stream){
		if(stream.buffer_data.pointer < 50){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Small: no compression... " << stream.buffer_data.pointer << std::endl;
			return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_data.pointer + 65536);
		//for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_data.data, stream.buffer_data.pointer, this->compression_level);
			if(ZSTD_isError(ret)){
				std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
				//exit(1);
			}
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << 8 << ": " << ret << '\t' << (ret == ZSTD_CONTENTSIZE_UNKNOWN) << '\t' << (ret == ZSTD_CONTENTSIZE_ERROR) << std::endl;
		//}

		const float fold = (float)stream.buffer_data.pointer/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Bad compression. Ignoring... " << stream.buffer_data.pointer << std::endl;

			/*
			// Debug
			if(stream.header.stride == 1){
				const BYTE* const d = reinterpret_cast<const BYTE* const>(stream.buffer_data.data);
				for(U32 i = 0; i < stream.n_entries; ++i)
					std::cerr << std::bitset<8>(d[i]) << ' ';
				std::cerr << std::endl;
			} else if(stream.header.stride == 2){
				const U16* const d = reinterpret_cast<const U16* const>(stream.buffer_data.data);
				for(U32 i = 0; i < stream.n_entries; ++i)
					std::cerr << std::bitset<16>(d[i]) << ' ';
				std::cerr << std::endl;
			}
			*/

			return true;
		}

		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Input: " << stream.buffer_data.pointer << " and output: " << ret << " -> " << (float)stream.buffer_data.pointer/ret << "-fold"  << std::endl;
		//std::cerr << "uniform: " << (int)stream.header.controller.uniform << '\t' << (int)stream.header.controller.type << std::endl;
		stream.header.uLength = stream.buffer_data.pointer;
		stream.header.cLength = ret;
		stream.header.controller.encoder = Core::ENCODE_ZSTD;
		memcpy(stream.buffer_data.data, this->buffer.data, ret);
		stream.buffer_data.pointer = ret;
		stream.header.n_extra = 1;
		stream.header.extra = new char[sizeof(BYTE)];
		stream.header.extra[0] = (BYTE)this->compression_level;

		return true;
	}

	const bool encodeStrides(stream_type& stream){
		if(stream.buffer_strides.pointer < 50){
			stream.header_stride.controller.encoder = Core::ENCODE_NONE;
			stream.header_stride.uLength = stream.buffer_strides.pointer;
			stream.header_stride.cLength = stream.buffer_strides.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Small stride: no compression... " << stream.buffer_data.pointer << std::endl;
			return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer + 65536);
		size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_strides.data, stream.buffer_strides.pointer, this->compression_level);
		if(ZSTD_isError(ret)){
			std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
		}

		const float fold = (float)stream.buffer_strides.pointer/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			stream.header_stride.controller.encoder = Core::ENCODE_NONE;
			stream.header_stride.uLength = stream.buffer_data.pointer;
			stream.header_stride.cLength = stream.buffer_data.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "STRIDE: Bad compression. Ignoring... " << stream.buffer_data.pointer << std::endl;
			return true;
		}

		std::cerr << Helpers::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << stream.buffer_strides.pointer << " and output: " << ret << " -> " << (float)stream.buffer_strides.pointer/ret << "-fold"  << std::endl;

		stream.header_stride.uLength = stream.buffer_strides.pointer;
		stream.header_stride.cLength = ret;
		stream.header_stride.controller.encoder = Core::ENCODE_ZSTD;
		memcpy(stream.buffer_strides.data, this->buffer.data, ret);
		stream.buffer_strides.pointer = ret;
		stream.header_stride.n_extra = 1;
		stream.header_stride.extra = new char[sizeof(BYTE)];
		stream.header_stride.extra[0] = (BYTE)this->compression_level;

		return true;
	}

	const bool encode(permutation_type& manager){
		this->buffer.reset();
		this->buffer.resize(manager.n_samples*sizeof(U32)+65536);

		// First
		const BYTE w = ceil(log2(manager.n_samples+1) / 8);
		if(w == 1){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (BYTE)manager[i];
		} else if(w == 2){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U16)manager[i];
			//manager.generateCRC();
			memset(manager.PPA.data, 0, this->buffer.pointer);
			manager.PPA.pointer = this->buffer.pointer;


		} else if(w == 3 || w == 4){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U32)manager[i];
		} else {
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U64)manager[i];
		}

		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer.data, this->buffer.pointer);
		manager.crc = crc;
		manager.u_length = this->buffer.pointer;
		bytePreprocessBits(&this->buffer.data[0], manager.PPA.pointer, &manager.PPA.data[0]);

		size_t ret = ZSTD_compress(this->buffer.data, this->buffer.pointer, manager.PPA.data, manager.PPA.pointer, this->compression_level);
		if(ZSTD_isError(ret)){
			std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
			//exit(1);
		}
		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "PPA in: " << manager.PPA.pointer << " and out: " << ret << std::endl;
		memcpy(manager.PPA.data, this->buffer.data, ret);
		manager.PPA.pointer = ret;
		manager.c_length = ret;

		return true;
	}

	const bool decode(stream_type& stream){ return false; }

private:
	int compression_level;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
