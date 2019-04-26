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

ZSTDCodec::~ZSTDCodec() {
	ZSTD_freeCCtx(this->compression_context_);
	ZSTD_freeDCtx(this->decompression_context_);
}

bool ZSTDCodec::Compress(const yon_buffer_t& src, yon_buffer_t& dst, const int compression_level) {
	dst.reset();
	dst.resize(src.size() + 65536);
	const size_t ret = ZSTD_compress(
							   dst.data(),
							   dst.capacity(),
							   src.data(),
							   src.size(),
							   compression_level);

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << src.size() << " and output: " << ret << " -> " << (float)src.size()/ret << "-fold"  << std::endl;

	if (ZSTD_isError(ret)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}
	dst.n_chars_ = ret;

	return true;
}

bool ZSTDCodec::Decompress(const yon_buffer_t& src, yon_buffer_t& dst) {
	const size_t ret = ZSTD_decompress(
							   dst.data(),
							   dst.capacity(),
							   src.data(),
							   src.size());

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << src.size() << " and output: " << ret << " -> " << (float)ret/src.size() << "-fold"  << std::endl;

	if (ZSTD_isError(ret)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	dst.n_chars_ = ret;

	return true;
}

bool ZSTDCodec::Compress(container_type& container) {
	if (container.header.n_entries == 0) {
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.data.n_chars_                         = 0;
		container.header.data_header.cLength            = 0;
		container.header.data_header.uLength            = 0;
	}

	if (container.header.data_header.controller.uniform || container.data_uncompressed.size() < 100) {
		container.data.resize(container.data_uncompressed.size() + 100);
		memcpy(container.data.data(), container.data_uncompressed.data(), container.data_uncompressed.size());
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.data.n_chars_                         = container.data_uncompressed.size();
		container.header.data_header.cLength            = container.data_uncompressed.size();
		container.header.data_header.uLength            = container.data_uncompressed.size();

		container.GenerateMd5();

		if (container.header.data_header.controller.mixedStride == true)
			return(this->CompressStrides(container));
		else return true;
	}

	this->buffer.reset();
	this->buffer.resize(container.data_uncompressed.size() + 65536);
	const size_t ret = ZSTD_compressCCtx(this->compression_context_,
							   this->buffer.data(),
							   this->buffer.capacity(),
							   container.data_uncompressed.data(),
							   container.data_uncompressed.size(),
							   this->compression_level_data);

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "Input: " << container.GetSizeUncompressed() << " and output: " << ret << " -> " << (float)container.GetSizeUncompressed()/ret << "-fold"  << std::endl;

	if (ZSTD_isError(ret)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	const float fold = (float)container.data_uncompressed.size() / ret;
	if (fold < MIN_COMPRESSION_FOLD) {
		container.data.resize(container.data_uncompressed.size() + 100);
		memcpy(container.data.data(), container.data_uncompressed.data(), container.data_uncompressed.size());
		container.header.data_header.controller.encoder = YON_ENCODE_NONE;
		container.data.n_chars_                         = container.data_uncompressed.size();
		container.header.data_header.cLength            = container.data_uncompressed.size();
		container.header.data_header.uLength            = container.data_uncompressed.size();

		container.GenerateMd5();

		if (container.header.data_header.controller.mixedStride == true)
			return(this->CompressStrides(container));
		else return true;
	}

	container.data.resize(ret + 65536);
	memcpy(container.data.data(), this->buffer.data(), ret);
	container.data.n_chars_                  = ret;
	container.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	container.header.data_header.cLength            = container.data.size();
	container.header.data_header.uLength            = container.data_uncompressed.size();

	container.GenerateMd5();

	if (container.header.data_header.controller.mixedStride == true)
		return(this->CompressStrides(container));
	else return true;
}

bool ZSTDCodec::CompressStrides(container_type& container) {
	if (container.header.stride_header.controller.uniform || container.strides_uncompressed.size() < 100) {
		container.strides.resize(container.strides_uncompressed.size() + 100);
		memcpy(container.strides.data(), container.strides_uncompressed.data(), container.strides_uncompressed.size());
		container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
		container.strides.n_chars_                        = container.strides_uncompressed.size();
		container.header.stride_header.cLength            = container.strides_uncompressed.size();
		container.header.stride_header.uLength            = container.strides_uncompressed.size();

		return true;
	}

	this->buffer.reset();
	this->buffer.resize(container.strides_uncompressed.size() + 65536);
	size_t ret = ZSTD_compressCCtx(this->compression_context_,
							   this->buffer.data(),
							   this->buffer.capacity(),
							   container.strides_uncompressed.data(),
							   container.strides_uncompressed.size(),
							   this->compression_level_data);

	if (ZSTD_isError(ret)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		return(false);
	}

	const float fold = (float)container.strides_uncompressed.size()/ret;
	if (fold < MIN_COMPRESSION_FOLD) {
		container.strides.resize(container.strides_uncompressed.size() + 100);
		memcpy(container.strides.data(), container.strides_uncompressed.data(), container.strides_uncompressed.size());
		container.header.stride_header.controller.encoder = YON_ENCODE_NONE;
		container.strides.n_chars_                        = container.strides_uncompressed.size();
		container.header.stride_header.cLength            = container.strides_uncompressed.size();
		container.header.stride_header.uLength            = container.strides_uncompressed.size();
		return true;
	}

	//std::cerr << utility::timestamp("LOG","COMPRESSION-STRIDE") << "Input: " << container.buffer_strides_uncompressed.n_chars << " and output: " << ret << " -> " << (float)container.buffer_strides_uncompressed.n_chars/ret << "-fold"  << std::endl;

	container.strides.resize(ret + 65536);
	memcpy(container.strides.data(), this->buffer.data(), ret);
	container.strides.n_chars_                        = ret;
	container.header.stride_header.cLength            = container.strides.size();
	container.header.stride_header.uLength            = container.strides_uncompressed.size();
	container.header.stride_header.controller.encoder = YON_ENCODE_ZSTD;

	return true;
}

bool ZSTDCodec::Compress(container_type& container, permutation_type& manager) {
	if (manager.n_s == 0)
		return true;

	container.data_uncompressed.reset();
	container.data_uncompressed.resize(manager.n_s*sizeof(uint32_t) + 65536);
	container.data_uncompressed << manager;

	this->buffer.reset();
	this->buffer.resize(manager.n_s*sizeof(uint32_t) + 65536);

	//const uint32_t in = manager.PPA.n_chars;
	const int p_ret = permuteIntBits(container.data_uncompressed.data(),
	                                 container.data_uncompressed.size(),
									 this->buffer.data());

	this->buffer.n_chars_ = p_ret;

	/*
	// DEBUG
	const int up_ret = unpermuteIntBits(this->buffer.data(),
										this->buffer.size(),
										manager.PPA.data());

	uint32_t crc2 = crc32(0, NULL, 0);
	crc2 = crc32(crc2, (Bytef*)manager.PPA.data(), up_ret);

	for (uint32_t i = 0; i < manager.n_samples; ++i)
		std::cerr << manager[i] << ' ';
	std::cerr << std::endl;

	std::cerr << in << "->" << p_ret << "->" << up_ret << std::endl;
	assert(manager.crc==crc2);
	*/


	container.data.reset();
	container.data.resize(p_ret + 65536);
	size_t ret = ZSTD_compress(container.data.data(),
	                           container.data.capacity(),
							   this->buffer.data(),
							   this->buffer.size(),
							   this->compression_level_data);


	if (ZSTD_isError(ret)) {
		std::cerr << "error zstd permute_ : " << ZSTD_getErrorCode(ret) << std::endl;
		std::cerr << ZSTD_getErrorName(ret) << std::endl;
		std::cerr << this->buffer.n_chars_ << '\t' << container.data.size() << std::endl;
		exit(1);
	}

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "PPA in: " << this->buffer.n_chars << " and out: " << ret << std::endl;

	container.data.n_chars_                         = ret;
	container.header.data_header.cLength            = container.data.size();
	container.header.data_header.uLength            = container.data_uncompressed.size();
	container.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	container.GenerateMd5();

	return true;
}

bool ZSTDCodec::Decompress(container_type& container) {
	if (container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE) {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if (container.header.data_header.controller.encoder != YON_ENCODE_ZSTD) {
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used in decompressor..." << std::endl;
		return false;
	}

	container.data_uncompressed.reset();
	container.data_uncompressed.resize(container.header.data_header.uLength + 65536);
	int ret = ZSTD_decompress(container.data_uncompressed.data(),
							  container.data_uncompressed.capacity(),
							  container.data.data(),
							  container.data.size());

	if (ZSTD_isError(ret)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
		std::cerr << ret << std::endl;
		return(false);
	}

	assert(ret >= 0);
	container.data_uncompressed.n_chars_ = ret;
	assert((uint32_t)ret == container.header.data_header.uLength);
	assert(container.CheckMd5(0));

	return true;
}

bool ZSTDCodec::DecompressStrides(container_type& container) {
	if (container.header.stride_header.controller.encryption != YON_ENCRYPTION_NONE) {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if (!container.header.data_header.controller.mixedStride)
		return false;

	if (container.header.stride_header.controller.encoder != YON_ENCODE_ZSTD) {
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used in decompressor..." << std::endl;
		return false;
	}

	container.strides_uncompressed.reset();
	container.strides_uncompressed.resize(container.header.stride_header.uLength + 65536);
	int ret_stride = ZSTD_decompress(container.strides_uncompressed.data(),
									 container.strides_uncompressed.capacity(),
									 container.strides.data(),
									 container.strides.size());

	if (ZSTD_isError(ret_stride)) {
		std::cerr << utility::timestamp("ERROR","ZSTD") << ZSTD_getErrorString(ZSTD_getErrorCode(ret_stride)) << std::endl;
		return(false);
	}

	assert(ret_stride >= 0);
	container.strides_uncompressed.n_chars_ = ret_stride;
	assert((uint32_t)ret_stride == container.header.stride_header.uLength);
	//std::cerr << "ENCODE_ZSTD | STRIDE | CRC check " << (container.checkCRC(0) ? "PASS" : "FAIL") << std::endl;
	assert(container.CheckMd5(1));

	return true;
}

bool ZSTDCodec::Decompress(container_type& container, permutation_type& manager) {
	this->buffer.reset();

	container.data_uncompressed.reset();
	container.data_uncompressed.resize(manager.n_s*sizeof(uint32_t) + 65536);
	this->buffer.resize(manager.n_s*sizeof(uint32_t) + 65536);
	size_t ret = ZSTD_decompress(this->buffer.data(),
								 this->buffer.capacity(),
								 container.data.data(),
								 container.data.size());

	if (ZSTD_isError(ret)) {
		std::cerr << "error zstd permute_ : " << ZSTD_getErrorCode(ret) << std::endl;
		std::cerr << ZSTD_getErrorName(ret) << std::endl;
		std::cerr << this->buffer.n_chars_ << '\t' << container.data_uncompressed.n_chars_ << std::endl;
		exit(1);
	}

	//std::cerr << utility::timestamp("LOG","COMPRESSION") << "PPA in: " << manager.PPA.size() << " and out: " << ret << std::endl;
	//manager.PPA.n_chars = ret;
	//manager.c_length    = ret;

	//container.buffer_data_uncompressed.resize(ret + 16536);
	const int up_ret = unpermuteIntBits(this->buffer.data(),
										ret,
										container.data_uncompressed.data());

		//std::cerr << "ret: " << up_ret << std::endl;
	//memcpy(manager.PPA.buffer, this->buffer.data(), up_ret);
	container.data_uncompressed.n_chars_ = up_ret;
	container.data_uncompressed >> manager;

	return true;
}

}
}
