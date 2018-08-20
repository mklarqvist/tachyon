#include "uncompressed_codec.h"

namespace tachyon{
namespace algorithm{

bool UncompressedCodec::Compress(container_type& container){
	container.buffer_data.resize(container.buffer_data_uncompressed.size() + 65536);
	memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
	container.header.data_header.controller.encoder = YON_ENCODE_NONE;
	container.buffer_data.n_chars_                  = container.buffer_data_uncompressed.size();
	container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
	return true;
}

bool UncompressedCodec::Decompress(container_type& container){
	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if(container.header.data_header.controller.encoder != YON_ENCODE_NONE){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
		return false;
	}
	container.buffer_data_uncompressed.resize(container.buffer_data.n_chars_ + 16536);
	memcpy(container.buffer_data_uncompressed.data(), container.buffer_data.data(), container.buffer_data.size());
	container.buffer_data_uncompressed.n_chars_ = container.buffer_data.size();
	assert(container.CheckMd5(0));
	return true;
}

bool UncompressedCodec::DecompressStrides(container_type& container){
	if(container.header.stride_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

	if(!container.header.data_header.controller.mixedStride){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Cannot decode strides. Stream has no strides..." << std::endl;
		return false;
	}

	if(container.header.stride_header.controller.encoder != YON_ENCODE_NONE){
		std::cerr << utility::timestamp("ERROR","ENCODER") << "Wrong codec used..." << std::endl;
		return false;
	}

	container.buffer_strides_uncompressed.resize(container.buffer_strides.n_chars_ + 16536);
	memcpy(container.buffer_strides_uncompressed.data(), container.buffer_strides.data(), container.buffer_strides.size());
	container.buffer_strides_uncompressed.n_chars_ = container.buffer_strides.n_chars_;
	assert(container.CheckMd5(1));
	return true;
}

}
}
