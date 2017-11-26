#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"
#include "zstd.h"
#include "common/zstd_errors.h"
#include "../../third_party/r16N.h"

namespace Tachyon{
namespace Compression{

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

	const bool encode(stream_type& stream){
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

		std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Input: " << stream.buffer_data.pointer << " and output: " << ret << " -> " << (float)stream.buffer_data.pointer/ret << "-fold"  << std::endl;
		stream.header.uLength = stream.buffer_data.pointer;
		stream.header.cLength = ret;
		stream.header.controller.encoder = Core::ENCODE_DEFLATE;
		//std::cerr << "DEFLATE STRIDE: " << stream.buffer_strides.pointer << '\t' << this->controller.buffer.size() << std::endl;

		memcpy(stream.buffer_data.data, this->buffer.data, ret);
		stream.buffer_data.pointer = ret;

		return true;
	}

	const int assessLevel(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_data.pointer + 65536);
		double prevRatio = 0;
		for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_data.data, stream.buffer_data.pointer, i);
			if(ZSTD_isError(ret)){
				std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
				//exit(1);
			}
			const double ratio = (double)stream.buffer_data.pointer/ret;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << i << ": " << stream.buffer_data.pointer << '\t' << ret << '\t' << ratio << std::endl;
			//if(ratio - prevRatio < 0.05) return(i);
			prevRatio = ratio;

		}
		return true;
	}

	const bool encodeStrides(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer + 65536);
		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "buffer size: " << this->buffer.capacity() << " and input: " << stream.buffer_data.pointer << std::endl;

		//for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_strides.data, stream.buffer_strides.pointer, this->compression_level);
			if(ZSTD_isError(ret)){
				std::cerr << "error: " << ZSTD_getErrorCode(ret) << std::endl;
				//exit(1);
			}
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << 8 << ": " << ret << '\t' << (ret == ZSTD_CONTENTSIZE_UNKNOWN) << '\t' << (ret == ZSTD_CONTENTSIZE_ERROR) << std::endl;
		//}

		std::cerr << Helpers::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << stream.buffer_strides.pointer << " and output: " << ret << " -> " << (float)stream.buffer_strides.pointer/ret << "-fold"  << std::endl;


		stream.header_stride.uLength = stream.buffer_strides.pointer;
		stream.header_stride.cLength = ret;
		stream.header_stride.controller.encoder = Core::ENCODE_DEFLATE;
		//std::cerr << "DEFLATE STRIDE: " << stream.buffer_strides.pointer << '\t' << this->controller.buffer.size() << std::endl;
		memcpy(stream.buffer_strides.data, this->buffer.data, ret);
		stream.buffer_strides.pointer = ret;

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
			memset(manager.PPA.data, 0, this->buffer.pointer);
			manager.PPA.pointer = this->buffer.pointer;
			const U32 block_size = this->buffer.pointer / w;
			std::cerr << "repacking to: " << block_size << " @ " << (int)w << std::endl;
			bytePreprocessBits(&this->buffer.data[0], manager.PPA.pointer, &manager.PPA.data[0]);

			//for(U32 i = 0; i < manager.PPA.pointer; ++i)
			//	std::cerr << std::bitset<8>(manager.PPA.data[i]);
			//std::cerr << std::endl;

		} else if(w == 3 || w == 4){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U32)manager[i];
		} else {
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U64)manager[i];
		}
		//memcpy(manager.PPA.data, this->buffer.data, this->buffer.pointer);
		//manager.PPA.pointer = this->buffer.pointer;
		//this->buffer.reset();

		//U32 out_size = this->buffer.capacity();
		//rans_compress_to_4x16((BYTE*)manager.PPA.data, manager.PPA.pointer, (BYTE*)this->buffer.data, &out_size, 1);
		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "PPA in: " << manager.PPA.pointer << " and out: " << out_size << std::endl;

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
