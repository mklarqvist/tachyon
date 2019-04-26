#include "uncompressed_codec.h"

namespace tachyon{
namespace algorithm{

bool UncompressedCodec::Compress(container_type& container) {
	container.data.resize(container.data_uncompressed.size() + 65536);
	memcpy(container.data.data(), container.data_uncompressed.data(), container.data_uncompressed.size());
	container.header.data_header.controller.encoder = YON_ENCODE_NONE;
	container.data.n_chars_                  = container.data_uncompressed.size();
	container.header.data_header.cLength            = container.data_uncompressed.size();
	return true;
}

bool UncompressedCodec::Decompress(container_type& container) {
	if (container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE) {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if (container.header.data_header.controller.encoder != YON_ENCODE_NONE) {
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
		return false;
	}
	container.data_uncompressed.resize(container.data.n_chars_ + 16536);
	memcpy(container.data_uncompressed.data(), container.data.data(), container.data.size());
	container.data_uncompressed.n_chars_ = container.data.size();
	assert(container.CheckMd5(0));
	return true;
}

bool UncompressedCodec::DecompressStrides(container_type& container) {
	if (container.header.stride_header.controller.encryption != YON_ENCRYPTION_NONE) {
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if (!container.header.data_header.controller.mixedStride) {
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Cannot decode strides. Stream has no strides..." << std::endl;
		return false;
	}

	if (container.header.stride_header.controller.encoder != YON_ENCODE_NONE) {
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
		return false;
	}

	container.strides_uncompressed.resize(container.strides.n_chars_ + 16536);
	memcpy(container.strides_uncompressed.data(), container.strides.data(), container.strides.size());
	container.strides_uncompressed.n_chars_ = container.strides.n_chars_;
	assert(container.CheckMd5(1));
	return true;
}

}
}
