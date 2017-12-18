#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"
#include "zstd.h"
#include "dictBuilder/zdict.h"
#include "common/zstd_errors.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"

namespace Tachyon{
namespace Compression{

// Lower bounds threshold in fold-change for compression
// to be kept
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
	virtual ~CompressionContainer(){ this->buffer.deleteAll(); }
	virtual const bool encode(permutation_type& manager) =0;
	virtual const bool encode(stream_type& stream) =0;
	virtual const bool encodeStrides(stream_type& stream) =0;
	virtual const bool decode(stream_type& stream) =0;
	virtual const bool decodeStrides(stream_type& stream) =0;

protected:
	buffer_type buffer;
};

class UncompressedCodec{
private:
	typedef UncompressedCodec self_type;
protected:
	typedef Core::StreamContainer stream_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::PermutationManager permutation_type;

public:
	UncompressedCodec(){}
	~UncompressedCodec(){}
	const bool encode(permutation_type& manager){ return false; }
	const bool encode(stream_type& stream){ return false; }
	const bool encodeStrides(stream_type& stream){ return false; }

	const bool decode(stream_type& stream){
		if(stream.header.controller.encoder != Core::ENCODE_NONE){
			std::cerr << "wrong codec" << std::endl;
			return false;
		}
		stream.buffer_data_uncompressed.resize(stream.buffer_data.pointer + 16536);
		memcpy(stream.buffer_data_uncompressed.data, stream.buffer_data.data, stream.buffer_data.pointer);
		stream.buffer_data_uncompressed.pointer = stream.buffer_data.pointer;
		assert(stream.checkCRC(0));
		return true;
	}

	const bool decodeStrides(stream_type& stream){
		if(!stream.header.controller.mixedStride){
			std::cerr << "has no mixed stride" << std::endl;
			return false;
		}

		if(stream.header_stride.controller.encoder != Core::ENCODE_NONE){
			std::cerr << "wrong codec" << std::endl;
			return false;
		}

		stream.buffer_strides_uncompressed.resize(stream.buffer_strides.pointer + 16536);
		memcpy(stream.buffer_strides_uncompressed.data, stream.buffer_strides.data, stream.buffer_strides.pointer);
		stream.buffer_strides_uncompressed.pointer = stream.buffer_strides.pointer;
		assert(stream.checkCRC(1));
		return true;
	}

protected:
	buffer_type buffer;
};

// GZIP of delfate codec
class ZSTDCodec : public CompressionContainer{
private:
	typedef ZSTDCodec self_type;
	typedef IO::TGZFController controller_type;

public:
	ZSTDCodec() : compression_level(0), cdict(nullptr){}
	~ZSTDCodec(){
		if(this->cdict != nullptr)
		ZSTD_freeCDict(this->cdict);
	}

	bool loadDictionary(const std::string file){
		std::ifstream stream(file, std::ios::binary | std::ios::ate);
		if(!stream.good()){
			std::cerr << "cannot open file: " << file << std::endl;
		}
		const size_t filesize = stream.tellg();
		stream.seekg(0);
		char* buffer = new char[filesize];
		stream.read(buffer, filesize);
		if(!stream.good()){
			std::cerr << "could not load data from: " << file << std::endl;
		}
		stream.close();

		this->cdict = ZSTD_createCDict(buffer, filesize, this->compression_level);
		return true;
	}

	void buildDictionary(buffer_type& buffer, buffer_type& dict_buffer, const U32& n_entries, const size_t* lengths){
		const size_t ret = ZDICT_trainFromBuffer(dict_buffer.data, dict_buffer.capacity(), buffer.data, lengths, n_entries);
		if(ZSTD_isError(ret)){
			std::cerr << "error zstd: " << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}
		dict_buffer.pointer = ret;
	}

	void setCompressionLevel(const int& c){ this->compression_level = c; }

	// not used
	const bool assess(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer + 65536);
		//for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_strides.data, stream.buffer_strides.pointer, this->compression_level);
			if(ZSTD_isError(ret)){
				std::cerr << "error zstd: " << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
				//exit(1);
			}
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << this->compression_level << ": " << stream.buffer_strides.pointer << '\t' << ret << '\t' << (double)stream.buffer_strides.pointer/ret << "-fold" << std::endl;
		//}

		//U32 out_size = this->buffer.capacity();
		//rans_compress_to_4x16((BYTE*)stream.buffer_strides.data, stream.buffer_strides.pointer, (BYTE*)this->buffer.data, &out_size, 0);
		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "RANS0\t" << stream.buffer_strides.pointer << "\t" << out_size << '\t' << (double)stream.buffer_strides.pointer/out_size << std::endl;

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
		stream.generateCRC();

		if(stream.header.controller.uniform || stream.buffer_data.pointer < 50){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			//stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Small: no compression... " << stream.buffer_data.pointer << std::endl;
			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_data.pointer + 65536);
		//for(S32 i = 1; i <= 22; ++i){
			size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_data.data, stream.buffer_data.pointer, this->compression_level);
			if(ZSTD_isError(ret)){
				std::cerr << "error zstd encode__: " << ZSTD_getErrorCode(ret) << std::endl;
				exit(1);
			}

			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "ZSTD@" << 8 << ": " << ret << '\t' << (ret == ZSTD_CONTENTSIZE_UNKNOWN) << '\t' << (ret == ZSTD_CONTENTSIZE_ERROR) << std::endl;
		//}

		const float fold = (float)stream.buffer_data.pointer/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			//stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Bad compression. Ignoring... " << stream.buffer_data.pointer << std::endl;

			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Input: " << stream.buffer_data.pointer << " and output: " << ret << " -> " << (float)stream.buffer_data.pointer/ret << "-fold"  << std::endl;
		//std::cerr << "uniform: " << (int)stream.header.controller.uniform << '\t' << (int)stream.header.controller.type << std::endl;
		//stream.header.uLength = stream.buffer_data.pointer;
		stream.header.cLength = ret;
		stream.header.controller.encoder = Core::ENCODE_ZSTD;
		memcpy(stream.buffer_data.data, this->buffer.data, ret);
		stream.buffer_data.pointer = ret;
		stream.header.n_extra = 1;
		stream.header.extra = new char[sizeof(BYTE)];
		stream.header.extra[0] = (BYTE)this->compression_level;

		if(stream.header.controller.mixedStride == true)
			return(this->encodeStrides(stream));
		else return true;
	}

	const bool encodeUsingDictionary(stream_type& stream){
		stream.generateCRC();

		if(stream.header.controller.uniform || stream.buffer_data.pointer < 50){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			//stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Small: no compression... " << stream.buffer_data.pointer << std::endl;
			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_data.pointer + 65536);
		//for(S32 i = 1; i <= 22; ++i){
		ZSTD_CCtx* cctx = ZSTD_createCCtx();
		size_t ret = ZSTD_compress_usingCDict(cctx,
											  this->buffer.data, this->buffer.capacity(),
											  stream.buffer_data.data, stream.buffer_data.pointer,
											  this->cdict);

		if(ZSTD_isError(ret)){
			std::cerr << "error zstd encode__: " << ZSTD_getErrorCode(ret) << std::endl;
			exit(1);
		}

		ZSTD_freeCCtx(cctx);

		const float fold = (float)stream.buffer_data.pointer/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			stream.header.controller.encoder = Core::ENCODE_NONE;
			//stream.header.uLength = stream.buffer_data.pointer;
			stream.header.cLength = stream.buffer_data.pointer;
			std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Bad compression. Ignoring... " << stream.buffer_data.pointer << std::endl;

			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Input: " << stream.buffer_data.pointer << " and output: " << ret << " -> " << (float)stream.buffer_data.pointer/ret << "-fold"  << std::endl;
		//std::cerr << "uniform: " << (int)stream.header.controller.uniform << '\t' << (int)stream.header.controller.type << std::endl;
		//stream.header.uLength = stream.buffer_data.pointer;
		stream.header.cLength = ret;
		stream.header.controller.encoder = Core::ENCODE_ZSTD;
		memcpy(stream.buffer_data.data, this->buffer.data, ret);
		stream.buffer_data.pointer = ret;
		stream.header.n_extra = 1;
		stream.header.extra = new char[sizeof(BYTE)];
		stream.header.extra[0] = (BYTE)this->compression_level;

		if(stream.header.controller.mixedStride == true)
			return(this->encodeStrides(stream));
		else return true;
	}


	const bool encodeStrides(stream_type& stream){
		if(stream.header_stride.controller.uniform || stream.buffer_strides.pointer < 50){
			stream.header_stride.controller.encoder = Core::ENCODE_NONE;
			//stream.header_stride.uLength = stream.buffer_strides.pointer;
			stream.header_stride.cLength = stream.buffer_strides.pointer;
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "Small stride: no compression... " << stream.buffer_data.pointer << std::endl;
			return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer + 65536);
		size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), stream.buffer_strides.data, stream.buffer_strides.pointer, this->compression_level);
		if(ZSTD_isError(ret)){
			std::cerr << "error zstd strides__: " << ZSTD_getErrorCode(ret) << std::endl;
			exit(1);
		}

		const float fold = (float)stream.buffer_strides.pointer/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			stream.header_stride.controller.encoder = Core::ENCODE_NONE;
			//stream.header_stride.uLength = stream.buffer_data.pointer;
			stream.header_stride.cLength = stream.buffer_data.pointer;
			//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "STRIDE: Bad compression. Ignoring... " << stream.buffer_data.pointer << std::endl;
			return true;
		}

		//std::cerr << Helpers::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << stream.buffer_strides.pointer << " and output: " << ret << " -> " << (float)stream.buffer_strides.pointer/ret << "-fold"  << std::endl;

		//stream.header_stride.uLength = stream.buffer_strides.pointer;
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
		this->buffer.resize(manager.n_samples*sizeof(U32) + 65536);

		// First
		const BYTE w = ceil(log2(manager.n_samples+1) / 8);
		if(w == 1){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (BYTE)manager[i];
		} else if(w == 2){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U16)manager[i];
		} else if(w == 3 || w == 4){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U32)manager[i];
		} else {
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U64)manager[i];
		}

		memset(manager.PPA.data, 0, this->buffer.pointer);
		manager.PPA.pointer = this->buffer.pointer;

		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer.data, this->buffer.pointer);
		manager.crc = crc;
		manager.u_length = this->buffer.pointer;
		bytePreprocessBits(&this->buffer.data[0], manager.PPA.pointer, &manager.PPA.data[0]);

		size_t ret = ZSTD_compress(this->buffer.data, this->buffer.capacity(), manager.PPA.data, manager.PPA.pointer, this->compression_level);
		if(ZSTD_isError(ret)){
			std::cerr << "error zstd permute_ : " << ZSTD_getErrorCode(ret) << std::endl;
			std::cerr << ZSTD_getErrorName(ret) << std::endl;
			std::cerr << this->buffer.pointer << '\t' << manager.PPA.pointer << std::endl;
			exit(1);
		}
		//std::cerr << Helpers::timestamp("LOG","COMPRESSION") << "PPA in: " << manager.PPA.pointer << " and out: " << ret << std::endl;
		memcpy(manager.PPA.data, this->buffer.data, ret);
		manager.PPA.pointer = ret;
		manager.c_length = ret;

		return true;
	}

	const bool decode(stream_type& stream){
		if(stream.header.controller.encoder != Core::ENCODE_ZSTD){
			std::cerr << "wrong codec" << std::endl;
			return false;
		}

		stream.buffer_data_uncompressed.resize(stream.header.uLength + 16536);
		int ret = ZSTD_decompress(stream.buffer_data_uncompressed.data,
								  stream.buffer_data_uncompressed.capacity(),
								  stream.buffer_data.data,
								  stream.buffer_data.pointer);

		assert(ret >= 0);
		stream.buffer_data_uncompressed.pointer = ret;
		assert((U32)ret == stream.header.uLength);
		assert(stream.checkCRC(0));
		//std::cerr << "ENCODE_ZSTD | CRC check " << (stream.checkCRC(0) ? "PASS" : "FAIL") << std::endl;

		/*
		if(stream.header.controller.mixedStride){
			this->decodeStrides(stream);
		}
		*/

		return true;
	}

	const bool decodeStrides(stream_type& stream){
		if(!stream.header.controller.mixedStride)
			return false;

		if(stream.header_stride.controller.encoder != Core::ENCODE_ZSTD){
			std::cerr << "wrong codec" << std::endl;
			return false;
		}


		if(stream.header_stride.controller.encoder != Core::ENCODE_ZSTD)
			return true;

		stream.buffer_strides_uncompressed.resize(stream.header_stride.uLength + 16536);
		int ret_stride = ZSTD_decompress(stream.buffer_strides_uncompressed.data,
								  stream.buffer_strides_uncompressed.capacity(),
								  stream.buffer_strides.data,
								  stream.buffer_strides.pointer);

		assert(ret_stride >= 0);
		stream.buffer_strides_uncompressed.pointer = ret_stride;
		assert((U32)ret_stride == stream.header_stride.uLength);
		//std::cerr << "ENCODE_ZSTD | STRIDE | CRC check " << (stream.checkCRC(0) ? "PASS" : "FAIL") << std::endl;
		assert(stream.checkCRC(1));

		return true;
	}

private:
	int compression_level;
	ZSTD_CDict* cdict;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
