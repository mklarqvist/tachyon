#include "zstd_codec.h"

#include "algorithm/digest/variant_digest_manager.h"

namespace tachyon{
namespace algorithm{

ZSTDCodec::ZSTDCodec() :
	compression_level_data(0),
	compression_level_strides(0),
	compression_context_(ZSTD_createCCtx()),
	decompression_context_(ZSTD_createDCtx())
{

}

ZSTDCodec::~ZSTDCodec(){
	ZSTD_freeCCtx(this->compression_context_);
	ZSTD_freeDCtx(this->decompression_context_);
}

bool ZSTDCodec::Compress(const io::BasicBuffer& src, io::BasicBuffer& dst, const int compression_level){
	dst.reset();
	dst.resize(src.size() + 65536);
	const size_t ret = ZSTD_compress(
							   dst.data(),
							   dst.capacity(),
							   src.data(),
							   src.size(),
							   compression_level);

	std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << src.size() << " and output: " << ret << " -> " << (float)src.size()/ret << "-fold"  << std::endl;

	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	return true;
}

bool ZSTDCodec::Compress(container_type& container){
	container.GenerateCRC();

	if(container.header.n_entries == 0){
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.buffer_data.n_chars                   = 0;
		container.header.data_header.cLength            = 0;
		container.header.data_header.uLength            = 0;
	}

	if(container.header.data_header.controller.uniform || container.buffer_data_uncompressed.size() < 100){
		memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
		container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
		container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

		if(container.header.data_header.controller.mixedStride == true)
			return(this->CompressStrides(container));
		else return true;
	}

	this->buffer.reset();
	this->buffer.resize(container.buffer_data_uncompressed.size() + 65536);
	const size_t ret = ZSTD_compressCCtx(this->compression_context_,
							   this->buffer.data(),
							   this->buffer.capacity(),
							   container.buffer_data_uncompressed.data(),
							   container.buffer_data_uncompressed.size(),
							   this->compression_level_data);

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << container.getSizeUncompressed() << " and output: " << ret << " -> " << (float)container.getSizeUncompressed()/ret << "-fold"  << std::endl;


	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	const float fold = (float)container.buffer_data_uncompressed.size() / ret;
	if(fold < MIN_COMPRESSION_FOLD){
		memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
		container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
		container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

		if(container.header.data_header.controller.mixedStride == true)
			return(this->CompressStrides(container));
		else return true;
	}

	container.buffer_data.resize(ret + 65536);
	memcpy(container.buffer_data.data(), this->buffer.data(), ret);
	container.buffer_data.n_chars                   = ret;
	container.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	container.header.data_header.cLength            = container.buffer_data.size();
	container.header.data_header.uLength            = container.buffer_data_uncompressed.size();

	if(container.header.data_header.controller.mixedStride == true)
		return(this->CompressStrides(container));
	else return true;
}

bool ZSTDCodec::CompressStrides(container_type& container){
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
	size_t ret = ZSTD_compressCCtx(this->compression_context_,
							   this->buffer.data(),
							   this->buffer.capacity(),
							   container.buffer_strides_uncompressed.data(),
							   container.buffer_strides_uncompressed.size(),
							   this->compression_level_data);

	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
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

	container.buffer_strides.resize(ret + 65536);
	memcpy(container.buffer_strides.data(), this->buffer.data(), ret);
	container.buffer_strides.n_chars                  = ret;
	container.header.stride_header.cLength            = container.buffer_strides.size();
	container.header.stride_header.uLength            = container.buffer_strides_uncompressed.size();
	container.header.stride_header.controller.encoder = YON_ENCODE_ZSTD;

	return true;
}

bool ZSTDCodec::Compress(permutation_type& manager){
	if(manager.PPA.size() == 0)
		return true;

	U32 n_samples = manager.n_samples;
	if(n_samples < 100){
		n_samples = 100;
		manager.PPA.resize(n_samples*sizeof(U32));
	}

	this->buffer.reset();
	this->buffer.resize(n_samples*sizeof(U32) + 65536);

	algorithm::VariantDigestManager::GenerateMd5(manager.PPA.data(), manager.PPA.size(), manager.header.data_header.crc);
	manager.header.data_header.uLength = manager.PPA.size();

	//const U32 in = manager.PPA.n_chars;
	const int p_ret = permuteIntBits(manager.PPA.data(),
									 manager.PPA.size(),
									 this->buffer.data());

	this->buffer.n_chars = p_ret;

	/*
	// DEBUG
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
	manager.header.data_header.uLength = this->buffer.n_chars;
	manager.header.data_header.cLength = ret;

	return true;
}

bool ZSTDCodec::Decompress(container_type& container){
	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if(container.header.data_header.controller.encoder != YON_ENCODE_ZSTD){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used in decompressor..." << std::endl;
		return false;
	}

	container.buffer_data_uncompressed.reset();
	container.buffer_data_uncompressed.resize(container.header.data_header.uLength + 65536);
	int ret = ZSTD_decompress(container.buffer_data_uncompressed.data(),
							  container.buffer_data_uncompressed.capacity(),
							  container.buffer_data.data(),
							  container.buffer_data.size());

	if(ZSTD_isError(ret)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		std::cerr << ret << std::endl;
		return(false);
	}

	assert(ret >= 0);
	container.buffer_data_uncompressed.n_chars = ret;
	assert((U32)ret == container.header.data_header.uLength);
	assert(container.CheckCRC(0));

	return true;
}

bool ZSTDCodec::DecompressStrides(container_type& container){
	if(container.header.stride_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if(!container.header.data_header.controller.mixedStride)
		return false;

	if(container.header.stride_header.controller.encoder != YON_ENCODE_ZSTD){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used in decompressor..." << std::endl;
		return false;
	}

	container.buffer_strides_uncompressed.reset();
	container.buffer_strides_uncompressed.resize(container.header.stride_header.uLength + 65536);
	int ret_stride = ZSTD_decompress(container.buffer_strides_uncompressed.data(),
									 container.buffer_strides_uncompressed.capacity(),
									 container.buffer_strides.data(),
									 container.buffer_strides.size());

	if(ZSTD_isError(ret_stride)){
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret_stride)) << std::endl;
		return(false);
	}

	assert(ret_stride >= 0);
	container.buffer_strides_uncompressed.n_chars = ret_stride;
	assert((U32)ret_stride == container.header.stride_header.uLength);
	//std::cerr << "ENCODE_ZSTD | STRIDE | CRC check " << (container.checkCRC(0) ? "PASS" : "FAIL") << std::endl;
	assert(container.CheckCRC(1));

	return true;
}

bool ZSTDCodec::Decompress(permutation_type& manager){
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

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "PPA in: " << manager.PPA.size() << " and out: " << ret << std::endl;
	//manager.PPA.n_chars = ret;
	//manager.c_length    = ret;

	manager.PPA.resize(ret + 16536);
	const int up_ret = unpermuteIntBits(this->buffer.data(),
										ret,
										manager.PPA.data());

	//std::cerr << "ret: " << up_ret << std::endl;
	//memcpy(manager.PPA.buffer, this->buffer.data(), up_ret);
	manager.PPA.n_chars = up_ret;
	return true;
}

}
}
