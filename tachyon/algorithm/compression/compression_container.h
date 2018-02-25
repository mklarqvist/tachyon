#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"
#include "zstd.h"
#include "zstd_errors.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"

namespace tachyon{
namespace algorithm{

/**< Lower bounds threshold in fold-change for compression to be kept */
#define MIN_COMPRESSION_FOLD 1.05

/**
 * Permute bits from a byte-stream of U32 into target
 * destinations such that bit X from a
 * byte [1,2,3,4,5,6,7,8] stream is permuted
 * to 1,1,1,1,1....8,8,8,8,8
 * @param data Input char* buffer
 * @param size Length of input data
 * @param destination Destination char* buffer of permuted data
 * @return TRUE if passing or FALSE otherwise
 */
inline const U32 permuteIntBits(const char* const  data,
                                        const U32  size,
                                             char* destination)
{
	if(size == 0) return 0;
	U32 internal_size = size + (32-size%32); // Balance bytes
	assert(internal_size % 32 == 0);

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	memset(dest, 0, internal_size); // Set all bytes to 0
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[32]; // Bucket pointers
	const U32 partition_size = internal_size / 32; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 32; ++i)
		target[31-i] = &dest[partition_size*i];

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i + 4 < internal_size; i+=4, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in U32
		// Update target T at byte position P with bit J at position K
		for(U32 j = 0; j < 8; ++j)
			target[j][p] |= ((d[i] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+8][p] |= ((d[i+1] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+16][p] |= ((d[i+2] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+24][p] |= ((d[i+3] & (1 << j)) >> j) << k;
	}

	return internal_size;
}

inline const U32 unpermuteIntBits(char* data,
                             const U32  size,
                                  char* destination)
{
	if(size == 0) return 0;
	//U32 internal_size = size + (32-size%32); // Balance bytes
	//assert(internal_size % 32 == 0);

	BYTE* temp = reinterpret_cast<BYTE*>(data); // Recast destination as U32
	U32* dest = reinterpret_cast<U32*>(destination); // Recast destination as U32
	const U32 n_entries = size / sizeof(U32);
	memset(destination, 0, size); // Set all bytes to 0
	//const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[32]; // Bucket pointers
	const U32 partition_size = size / 32; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 32; ++i)
		target[31-i] = &temp[partition_size*i];

	/*
	for(U32 i = 0; i < size; ++i){
		std::cerr << (int)data[i] << ' ';
	}
	std::cerr << std::endl;
	*/

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i < n_entries; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(U32 j = 0; j < 32; ++j){
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
		}
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(U32) << std::endl;
	return size;
}

class CompressionContainer{
private:
	typedef CompressionContainer self_type;
protected:
	typedef containers::DataContainer stream_type;
	typedef io::BasicBuffer buffer_type;
	typedef algorithm::PermutationManager permutation_type;

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
	typedef UncompressedCodec             self_type;

protected:
	typedef containers::DataContainer     stream_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	UncompressedCodec(){}
	~UncompressedCodec(){ this->buffer.deleteAll(); }
	inline const bool encode(permutation_type& manager){ return true; }
	inline const bool encode(stream_type& stream){ return true; }
	inline const bool encodeStrides(stream_type& stream){ return true; }

	const bool decode(stream_type& stream){
		if(stream.header.controller.encoder != core::YON_ENCODE_NONE){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}
		stream.buffer_data_uncompressed.resize(stream.buffer_data.n_chars + 16536);
		memcpy(stream.buffer_data_uncompressed.buffer, stream.buffer_data.buffer, stream.buffer_data.n_chars);
		stream.buffer_data_uncompressed.n_chars = stream.buffer_data.n_chars;
		assert(stream.checkCRC(0));
		return true;
	}

	const bool decodeStrides(stream_type& stream){
		if(!stream.header.controller.mixedStride){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Cannot decode strides. Stream has no strides..." << std::endl;
			return false;
		}

		if(stream.header_stride.controller.encoder != core::YON_ENCODE_NONE){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}

		stream.buffer_strides_uncompressed.resize(stream.buffer_strides.n_chars + 16536);
		memcpy(stream.buffer_strides_uncompressed.buffer, stream.buffer_strides.buffer, stream.buffer_strides.n_chars);
		stream.buffer_strides_uncompressed.n_chars = stream.buffer_strides.n_chars;
		assert(stream.checkCRC(1));
		return true;
	}

protected:
	buffer_type buffer;
};

// ZSTD codec
class ZSTDCodec : public CompressionContainer{
private:
	typedef ZSTDCodec          self_type;
	typedef io::TGZFController controller_type;

public:
	ZSTDCodec() : compression_level_data(0), compression_level_strides(0){}
	~ZSTDCodec(){}

	inline void setCompressionLevel(const int& c){ this->compression_level_data = c; this->compression_level_strides = c; }
	inline void setCompressionLevelData(const int& c){ this->compression_level_data = c; }
	inline void setCompressionLevelStrides(const int& c){ this->compression_level_strides = c; }

	/**
	 *
	 * @param stream
	 * @return
	 */
	const bool encode(stream_type& stream){
		stream.generateCRC();

		if(stream.header.controller.uniform || stream.buffer_data_uncompressed.size() < 50){
			memcpy(stream.buffer_data.data(), stream.buffer_data_uncompressed.data(), stream.buffer_data_uncompressed.size());
			stream.header.controller.encoder = core::YON_ENCODE_NONE;
			stream.buffer_data.n_chars       = stream.buffer_data_uncompressed.size();
			stream.header.cLength            = stream.buffer_data_uncompressed.size();

			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_data_uncompressed.size() + 65536);
		size_t ret = ZSTD_compress(this->buffer.data(),
								   this->buffer.capacity(),
								   stream.buffer_data_uncompressed.data(),
								   stream.buffer_data_uncompressed.size(),
								   this->compression_level_data);

		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		const float fold = (float)stream.buffer_data_uncompressed.size() / ret;
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(stream.buffer_data.data(), stream.buffer_data_uncompressed.data(), stream.buffer_data_uncompressed.size());
			stream.header.controller.encoder = core::YON_ENCODE_NONE;
			stream.buffer_data.n_chars       = stream.buffer_data_uncompressed.size();
			stream.header.cLength            = stream.buffer_data_uncompressed.size();

			if(stream.header.controller.mixedStride == true)
				return(this->encodeStrides(stream));
			else return true;
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << stream.buffer_data.n_chars << " and output: " << ret << " -> " << (float)stream.buffer_data.n_chars/ret << "-fold"  << std::endl;

		memcpy(stream.buffer_data.data(), this->buffer.data(), ret);
		stream.header.cLength            = ret;
		stream.header.controller.encoder = core::YON_ENCODE_ZSTD;
		stream.buffer_data.n_chars       = ret;
		stream.header.n_extra            = 1;
		stream.header.extra              = new char[sizeof(BYTE)];
		stream.header.extra[0]           = (BYTE)this->compression_level_data;

		if(stream.header.controller.mixedStride == true)
			return(this->encodeStrides(stream));
		else return true;
	}

	/**
	 *
	 * @param stream
	 * @return
	 */
	const bool encodeStrides(stream_type& stream){
		if(stream.header_stride.controller.uniform || stream.buffer_strides_uncompressed.size() < 50){
			memcpy(stream.buffer_strides.data(), stream.buffer_strides_uncompressed.data(), stream.buffer_strides_uncompressed.size());
			stream.header_stride.controller.encoder = core::YON_ENCODE_NONE;
			stream.buffer_strides.n_chars           = stream.buffer_strides_uncompressed.size();
			stream.header_stride.cLength            = stream.buffer_strides_uncompressed.size();

			return true;
		}

		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides_uncompressed.size() + 65536);
		size_t ret = ZSTD_compress(this->buffer.data(),
                                   this->buffer.capacity(),
								   stream.buffer_strides_uncompressed.data(),
								   stream.buffer_strides_uncompressed.size(),
								   this->compression_level_data);

		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		const float fold = (float)stream.buffer_strides_uncompressed.size()/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(stream.buffer_strides.data(), stream.buffer_strides_uncompressed.data(), stream.buffer_strides_uncompressed.size());
			stream.header_stride.controller.encoder = core::YON_ENCODE_NONE;
			stream.buffer_strides.n_chars           = stream.buffer_strides_uncompressed.size();
			stream.header_stride.cLength            = stream.buffer_strides_uncompressed.size();
			return true;
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << stream.buffer_strides_uncompressed.n_chars << " and output: " << ret << " -> " << (float)stream.buffer_strides_uncompressed.n_chars/ret << "-fold"  << std::endl;

		memcpy(stream.buffer_strides.data(), this->buffer.data(), ret);
		stream.header_stride.cLength            = ret;
		stream.header_stride.controller.encoder = core::YON_ENCODE_ZSTD;
		stream.buffer_strides.n_chars           = ret;
		stream.header_stride.n_extra            = 1;
		stream.header_stride.extra              = new char[sizeof(BYTE)];
		stream.header_stride.extra[0]           = (BYTE)this->compression_level_data;

		return true;
	}

	/**
	 *
	 * @param manager
	 * @return
	 */
	const bool encode(permutation_type& manager){
		if(manager.PPA.size() == 0)
			return true;

		U32 n_samples = manager.n_samples;
		if(n_samples < 100){
			n_samples = 100;
			manager.PPA.resize(n_samples*sizeof(U32));
		}

		this->buffer.reset();
		this->buffer.resize(n_samples*sizeof(U32) + 65536);

		U32 crc          = crc32(0, NULL, 0);
		manager.crc      = crc32(crc, (Bytef*)manager.PPA.data(), manager.PPA.size());
		manager.u_length = manager.PPA.size();

		//const U32 in = manager.PPA.n_chars;
		const int p_ret = permuteIntBits(manager.PPA.data(),
				                         manager.PPA.size(),
										 this->buffer.data());

		this->buffer.n_chars = p_ret;

		/*
		const int up_ret = unpermuteIntBits(this->buffer.data(),
											this->buffer.size(),
										    manager.PPA.data());

		U32 crc2 = crc32(0, NULL, 0);
		crc2 = crc32(crc2, (Bytef*)manager.PPA.data(), up_ret);

		for(U32 i = 0; i < manager.n_samples; ++i)
			std::cerr << manager[i] << ' ';
		std::cerr << std::endl;

		std::cerr << in << "->" << p_ret << "->" << up_ret << std::endl;
		assert(manager.crc==crc2);
		*/


		size_t ret = ZSTD_compress(manager.PPA.data(),
								   manager.PPA.capacity(),
								   this->buffer.data(),
								   this->buffer.size(),
								   this->compression_level_data);


		if(ZSTD_isError(ret)){
			std::cerr << "error zstd permute_ : " << ZSTD_getErrorCode(ret) << std::endl;
			std::cerr << ZSTD_getErrorName(ret) << std::endl;
			std::cerr << this->buffer.n_chars << '\t' << manager.PPA.n_chars << std::endl;
			exit(1);
		}
		//std::cerr << utility::timestamp("LOG","COMPRESSION") << "PPA in: " << this->buffer.n_chars << " and out: " << ret << std::endl;
		manager.PPA.n_chars = ret;
		manager.c_length    = ret;

		return true;
	}

	const bool decode(stream_type& stream){
		if(stream.header.controller.encoder != core::YON_ENCODE_ZSTD){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}

		stream.buffer_data_uncompressed.resize(stream.header.uLength + 16536);
		int ret = ZSTD_decompress(stream.buffer_data_uncompressed.data(),
								  stream.buffer_data_uncompressed.capacity(),
								  stream.buffer_data.data(),
								  stream.buffer_data.size());

		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		assert(ret >= 0);
		stream.buffer_data_uncompressed.n_chars = ret;
		assert((U32)ret == stream.header.uLength);
		assert(stream.checkCRC(0));

		return true;
	}

	const bool decodeStrides(stream_type& stream){
		if(!stream.header.controller.mixedStride)
			return false;

		if(stream.header_stride.controller.encoder != core::YON_ENCODE_ZSTD){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}


		if(stream.header_stride.controller.encoder != core::YON_ENCODE_ZSTD)
			return true;

		stream.buffer_strides_uncompressed.resize(stream.header_stride.uLength + 16536);
		int ret_stride = ZSTD_decompress(stream.buffer_strides_uncompressed.data(),
										 stream.buffer_strides_uncompressed.capacity(),
										 stream.buffer_strides.data(),
										 stream.buffer_strides.size());

		assert(ret_stride >= 0);
		stream.buffer_strides_uncompressed.n_chars = ret_stride;
		assert((U32)ret_stride == stream.header_stride.uLength);
		//std::cerr << "ENCODE_ZSTD | STRIDE | CRC check " << (stream.checkCRC(0) ? "PASS" : "FAIL") << std::endl;
		assert(stream.checkCRC(1));

		return true;
	}

	const bool decode(permutation_type& manager){
		this->buffer.reset();
		U32 n_samples = manager.n_samples;
		if(n_samples < 100){
			n_samples = 100;
			manager.PPA.resize(n_samples*sizeof(U32));
		}
		this->buffer.resize(n_samples*sizeof(U32) + 65536);
		size_t ret = ZSTD_decompress(this->buffer.data(),
			    					 this->buffer.capacity(),
									 manager.PPA.data(),
								     manager.PPA.size());

		if(ZSTD_isError(ret)){
			std::cerr << "error zstd permute_ : " << ZSTD_getErrorCode(ret) << std::endl;
			std::cerr << ZSTD_getErrorName(ret) << std::endl;
			std::cerr << this->buffer.n_chars << '\t' << manager.PPA.n_chars << std::endl;
			exit(1);
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION") << "PPA in: " << this->buffer.n_chars << " and out: " << ret << std::endl;
		//manager.PPA.n_chars = ret;
		//manager.c_length    = ret;

		manager.PPA.resize(ret + 16536);
		const int up_ret = unpermuteIntBits(this->buffer.data(),
											ret,
											manager.PPA.data());

		//memcpy(manager.PPA.buffer, this->buffer.data(), up_ret);
		manager.PPA.n_chars = up_ret;
		return true;
	}

private:
	int compression_level_data;
	int compression_level_strides;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
