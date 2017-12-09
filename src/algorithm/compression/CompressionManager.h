#ifndef ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_

#include "CompressionContainer.h"
#include "../../core/StreamContainer.h"

namespace Tachyon{
namespace Compression{

class CompressionManager{
private:
	typedef CompressionManager self_type;
	typedef UncompressedCodec no_codec_type;
	typedef ZSTDCodec zstd_codec_type;
	typedef Core::StreamContainer container_type;

public:
	CompressionManager(){}
	~CompressionManager(){}

	bool decompress(container_type& container){
		if(container.header.controller.encoder == Core::ENCODE_ZSTD){
			if(!this->zstd_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else if(container.header.controller.encoder == Core::ENCODE_NONE){
			if(!this->no_codec.decode(container)){ std::cerr << "failed" << std::endl; return false; }
		} else {
			std::cerr << "ILLEGAL CODEC" << std::endl;
			return false;
		}

		if(container.header.controller.mixedStride){
			if(container.header_stride.controller.encoder == Core::ENCODE_ZSTD){
				this->zstd_codec.decodeStrides(container);
			} else if (container.header_stride.controller.encoder == Core::ENCODE_NONE){
				this->no_codec.decodeStrides(container);
			}
		}
		return true;
	}

private:
	no_codec_type no_codec;
	zstd_codec_type zstd_codec;
};

}
}



#endif /* ALGORITHM_COMPRESSION_COMPRESSIONMANAGER_H_ */
