#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"
#include "zstd.h"
#include "zstd_errors.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"
#include "zpaq_wrapper.h"

namespace tachyon{
namespace algorithm{

/**< Lower bounds threshold in fold-change for compression to be kept */
#define MIN_COMPRESSION_FOLD 1.1

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

inline const U32 permuteByteBits(const char* const  data,
                                        const U32  size,
                                             char* destination)
{
	if(size == 0) return 0;
	U32 internal_size = size + (8 - size % 8); // Balance bytes
	assert(internal_size % 8 == 0);

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	memset(dest, 0, internal_size); // Set all bytes to 0
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[8]; // Bucket pointers
	const U32 partition_size = internal_size / 8; // Partition size
	assert(partition_size != 0);

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &dest[partition_size*i];

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i < internal_size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in U32
		// Update bucket B at byte position P with bit J at position K
		for(U32 bucket = 0; bucket < 8; ++bucket)
			target[bucket][p] |= ((d[i] & (1 << bucket)) >> bucket) << k;
	}

	return internal_size;
}

inline const U32 unpermuteByteBits(char* data,
                             const U32  size,
                                  char* destination)
{
	if(size == 0) return 0;
	//U32 internal_size = size + (32-size%32); // Balance bytes
	//assert(internal_size % 32 == 0);

	BYTE* temp = reinterpret_cast<BYTE*>(data); // Recast destination as U32
	BYTE* dest = reinterpret_cast<BYTE*>(destination); // Recast destination as U32
	memset(destination, 0, size); // Set all bytes to 0
	//const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[8]; // Bucket pointers
	const U32 partition_size = size / 8; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &temp[partition_size*i];

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
	for(U32 i = 0; i < size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(U32 j = 0; j < 8; ++j){
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
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	CompressionContainer(){}
	virtual ~CompressionContainer(){ }
	virtual const bool compress(permutation_type& manager) =0;
	virtual const bool compress(container_type& container) =0;
	virtual const bool compressStrides(container_type& container) =0;
	virtual const bool decompress(container_type& container) =0;
	virtual const bool decompressStrides(container_type& container) =0;

protected:
	buffer_type buffer;
};

class UncompressedCodec{
private:
	typedef UncompressedCodec             self_type;

protected:
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	UncompressedCodec(){}
	~UncompressedCodec(){ }
	inline const bool compress(permutation_type& manager){ return true; }
	inline const bool compress(container_type& container){
		container.buffer_data.resize(container.buffer_data_uncompressed.size() + 65536);
		memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.buffer_data.n_chars       = container.buffer_data_uncompressed.size();
		container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
		return true;
	}
	inline const bool compressStrides(container_type& container){ return true; }

	const bool decompress(container_type& container){
		if(container.header.data_header.controller.encoder != YON_ENCODE_NONE){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}
		container.buffer_data_uncompressed.resize(container.buffer_data.n_chars + 16536);
		memcpy(container.buffer_data_uncompressed.buffer, container.buffer_data.buffer, container.buffer_data.n_chars);
		container.buffer_data_uncompressed.n_chars = container.buffer_data.n_chars;
		assert(container.checkCRC(0));
		return true;
	}

	const bool decompressStrides(container_type& container){
		if(!container.header.data_header.controller.mixedStride){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Cannot decode strides. Stream has no strides..." << std::endl;
			return false;
		}

		if(container.header.stride_header.controller.encoder != YON_ENCODE_NONE){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}

		container.buffer_strides_uncompressed.resize(container.buffer_strides.n_chars + 16536);
		memcpy(container.buffer_strides_uncompressed.buffer, container.buffer_strides.buffer, container.buffer_strides.n_chars);
		container.buffer_strides_uncompressed.n_chars = container.buffer_strides.n_chars;
		assert(container.checkCRC(1));
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
	const bool compress(container_type& container){
		container.generateCRC();

		if(container.header.data_header.controller.uniform || container.buffer_data_uncompressed.size() < 100){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
			container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

			if(container.header.data_header.controller.mixedStride == true)
				return(this->compressStrides(container));
			else return true;
		}

		this->buffer.reset();
		this->buffer.resize(container.buffer_data_uncompressed.size() + 65536);
		size_t ret = ZSTD_compress(this->buffer.data(),
								   this->buffer.capacity(),
								   container.buffer_data_uncompressed.data(),
								   container.buffer_data_uncompressed.size(),
								   this->compression_level_data);

		//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << container.getSizeUncompressed() << " and output: " << ret << " -> " << (float)container.getSizeUncompressed()/ret << "-fold"  << std::endl;


		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		const float fold = (float)container.buffer_data_uncompressed.size() / ret;
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
			container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

			if(container.header.data_header.controller.mixedStride == true)
				return(this->compressStrides(container));
			else return true;
		}

		container.buffer_data.resize(ret + 65536);
		memcpy(container.buffer_data.data(), this->buffer.data(), ret);
		container.header.data_header.cLength            = ret;
		container.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
		container.buffer_data.n_chars                   = ret;
		container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

		if(container.header.data_header.controller.mixedStride == true)
			return(this->compressStrides(container));
		else return true;
	}

	/**
	 *
	 * @param stream
	 * @return
	 */
	const bool compressStrides(container_type& container){
		if(container.header.stride_header.controller.uniform || container.buffer_strides_uncompressed.size() < 100){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars                  = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();
			container.header.stride_header.uLength            = container.buffer_strides_uncompressed.size();

			return true;
		}

		this->buffer.reset();
		this->buffer.resize(container.buffer_strides_uncompressed.size() + 65536);
		size_t ret = ZSTD_compress(this->buffer.data(),
                                   this->buffer.capacity(),
								   container.buffer_strides_uncompressed.data(),
								   container.buffer_strides_uncompressed.size(),
								   this->compression_level_data);

		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		const float fold = (float)container.buffer_strides_uncompressed.size()/ret;
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars                  = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();
			container.header.stride_header.uLength            = container.buffer_strides_uncompressed.size();
			return true;
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << container.buffer_strides_uncompressed.n_chars << " and output: " << ret << " -> " << (float)container.buffer_strides_uncompressed.n_chars/ret << "-fold"  << std::endl;

		container.buffer_data.resize(ret + 65536);
		memcpy(container.buffer_strides.data(), this->buffer.data(), ret);
		container.header.stride_header.cLength            = ret;
		container.header.stride_header.controller.encoder = YON_ENCODE_ZSTD;
		container.buffer_strides.n_chars                  = ret;

		return true;
	}

	/**
	 *
	 * @param manager
	 * @return
	 */
	const bool compress(permutation_type& manager){
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

	const bool decompress(container_type& container){
		if(container.header.data_header.controller.encoder != YON_ENCODE_ZSTD){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}
		container.buffer_data_uncompressed.resize(container.header.data_header.uLength + 16536);
		int ret = ZSTD_decompress(container.buffer_data_uncompressed.data(),
								  container.buffer_data_uncompressed.capacity(),
								  container.buffer_data.data(),
								  container.buffer_data.size());

		if(ZSTD_isError(ret)){
			std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			exit(1);
		}

		assert(ret >= 0);
		container.buffer_data_uncompressed.n_chars = ret;
		assert((U32)ret == container.header.data_header.uLength);
		assert(container.checkCRC(0));

		return true;
	}

	const bool decompressStrides(container_type& container){
		if(!container.header.data_header.controller.mixedStride)
			return false;

		if(container.header.stride_header.controller.encoder != YON_ENCODE_ZSTD){
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
			return false;
		}


		if(container.header.stride_header.controller.encoder != YON_ENCODE_ZSTD)
			return true;

		container.buffer_strides_uncompressed.resize(container.header.stride_header.uLength + 16536);
		int ret_stride = ZSTD_decompress(container.buffer_strides_uncompressed.data(),
										 container.buffer_strides_uncompressed.capacity(),
										 container.buffer_strides.data(),
										 container.buffer_strides.size());

		assert(ret_stride >= 0);
		container.buffer_strides_uncompressed.n_chars = ret_stride;
		assert((U32)ret_stride == container.header.stride_header.uLength);
		//std::cerr << "ENCODE_ZSTD | STRIDE | CRC check " << (container.checkCRC(0) ? "PASS" : "FAIL") << std::endl;
		assert(container.checkCRC(1));

		return true;
	}

	const bool decompress(permutation_type& manager){
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

// ZPAQ
class ZPAQContainer : public CompressionContainer{
private:
	typedef ZPAQContainer self_type;

public:
	ZPAQContainer() :
		compression_level_data(3),
		compression_level_strides(3)
	{
	}

	virtual ~ZPAQContainer(){ }

	const bool compress(permutation_type& manager){ return false; }

	const bool compress(container_type& container, const std::string& command, const bool compress_strides = true){
		container.generateCRC();

		if(container.header.data_header.controller.uniform || container.buffer_data_uncompressed.size() < 100){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars       = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.controller.mixedStride == true)
					return(this->compressStrides(container, command));
				else return true;
			} else return true;
		}


		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, &command[0], "1", NULL, false);

		const float fold = (float)container.buffer_data_uncompressed.size() / out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars       = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
			if(compress_strides){
				if(container.header.data_header.controller.mixedStride == true)
					return(this->compressStrides(container, command));
				else return true;
			} else return true;
		}

		container.header.data_header.cLength            = out.buffer.size();
		container.header.data_header.controller.encoder = YON_ENCODE_ZPAQ;
		if(compress_strides){
			if(container.header.data_header.controller.mixedStride == true)
				return(this->compressStrides(container, command));
			else return true;
		} else return true;
	}

	const bool compress(container_type& container){
		return(this->compress(container, true));
	}

	const bool compress(container_type& container, const bool compress_strides){
		container.generateCRC();

		if(container.header.data_header.controller.uniform || container.buffer_data_uncompressed.size() < 100){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars       = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.controller.mixedStride == true)
					return(this->compressStrides(container));
				else return true;
			} else return true;
			//return true;
		}


		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, "x0.3ci1", "1", NULL, false);

		const float fold = (float)container.buffer_data_uncompressed.size() / out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars       = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.controller.mixedStride == true)
					return(this->compressStrides(container));
				else return true;
			} else return true;
		}

		container.header.data_header.cLength            = out.buffer.size();
		container.header.data_header.controller.encoder = YON_ENCODE_ZPAQ;
		if(compress_strides){
			if(container.header.data_header.controller.mixedStride == true)
				return(this->compressStrides(container));
			else return true;
		} else return true;
	}

	const bool compressStrides(container_type& container, const std::string& command){
		if(container.header.stride_header.controller.uniform || container.buffer_strides_uncompressed.size() < 100){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars                  = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();

			return true;
		}

		container.buffer_data.reset();
		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, &command[0], "2", NULL, false);

		const float fold = (float)container.buffer_strides_uncompressed.size()/out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars                  = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();
			return true;
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << container.buffer_strides_uncompressed.n_chars << " and output: " << ret << " -> " << (float)container.buffer_strides_uncompressed.n_chars/ret << "-fold"  << std::endl;

		container.header.stride_header.cLength            = out.buffer.size();
		container.header.stride_header.controller.encoder = YON_ENCODE_ZPAQ;

		return true;
	}

	const bool compressStrides(container_type& container){
		if(container.header.stride_header.controller.uniform || container.buffer_strides_uncompressed.size() < 100){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars           = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();

			return true;
		}

		container.buffer_data.reset();
		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, "x0.3ci1", "2", NULL, false);

		const float fold = (float)container.buffer_strides_uncompressed.size()/out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_strides.data(), container.buffer_strides_uncompressed.data(), container.buffer_strides_uncompressed.size());
			container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_strides.n_chars           = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();
			return true;
		}

		//std::cerr << utility::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << container.buffer_strides_uncompressed.n_chars << " and output: " << ret << " -> " << (float)container.buffer_strides_uncompressed.n_chars/ret << "-fold"  << std::endl;

		container.header.stride_header.cLength            = out.buffer.size();
		container.header.stride_header.controller.encoder = YON_ENCODE_ZPAQ;

		return true;
	}
	const bool decompress(container_type& container){
		if(container.header.data_header.controller.encoder != YON_ENCODE_ZPAQ){
			return true;
		}

		container.buffer_data_uncompressed.reset();
		container.buffer_data_uncompressed.resize(container.header.data_header.uLength + 65536);
		ZpaqWrapperIn in(container.buffer_data);
		ZpaqWrapperOut out(container.buffer_data_uncompressed);
		libzpaq::decompress(&in, &out);
		assert(out.buffer.size() == container.header.data_header.uLength);
		assert(container.checkCRC(0));

		return true;
	}
	const bool decompressStrides(container_type& container){
		if(container.header.stride_header.controller.encoder != YON_ENCODE_ZPAQ){
			return true;
		}
		std::cerr << "at strides" << std::endl;

		container.buffer_strides_uncompressed.reset();
		container.buffer_strides_uncompressed.resize(container.header.stride_header.uLength + 65536);
		ZpaqWrapperIn in(container.buffer_strides);
		ZpaqWrapperOut out(container.buffer_strides_uncompressed);
		libzpaq::decompress(&in, &out);
		assert(out.buffer.size() == container.header.stride_header.uLength);
		assert(container.checkCRC(1));

		return true;
	}

protected:
	int compression_level_data;
	int compression_level_strides;
	std::string compression_level_data_string;
	std::string compression_level_strides_string;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
