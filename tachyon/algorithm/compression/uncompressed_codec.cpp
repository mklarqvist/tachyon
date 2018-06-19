#include "uncompressed_codec.h"

namespace tachyon{
namespace algorithm{

const bool UncompressedCodec::compress(container_type& container){
	container.buffer_data.resize(container.buffer_data_uncompressed.size() + 65536);
	memcpy(container.buffer_data.data(), container.buffer_data_uncompressed.data(), container.buffer_data_uncompressed.size());
	container.header.data_header.controller.encoder = YON_ENCODE_NONE;
	container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
	container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
	return true;
}

const bool UncompressedCodec::decompress(container_type& container){
	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_NONE){
		std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Data is encrypted. Provide a valid keychain and decrypt before proceeding..." << std::endl;
		return false;
	}

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

const bool UncompressedCodec::decompressStrides(container_type& container){
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

	container.buffer_strides_uncompressed.resize(container.buffer_strides.n_chars + 16536);
	memcpy(container.buffer_strides_uncompressed.buffer, container.buffer_strides.buffer, container.buffer_strides.n_chars);
	container.buffer_strides_uncompressed.n_chars = container.buffer_strides.n_chars;
	assert(container.checkCRC(1));
	return true;
}

}
}
