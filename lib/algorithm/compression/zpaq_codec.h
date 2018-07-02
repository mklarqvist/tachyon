#ifndef ALGORITHM_COMPRESSION_ZPAQ_CODEC_H_
#define ALGORITHM_COMPRESSION_ZPAQ_CODEC_H_

#include "compression_container.h"
#include "zpaq_wrapper.h"

namespace tachyon{
namespace algorithm{

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
			memcpy(container.buffer_data.data(),
                   container.buffer_data_uncompressed.data(),
                   container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.hasMixedStride())
					return(this->compressStrides(container, command));
				else return true;
			} else return true;
		}

		container.buffer_data.reset();
		container.buffer_data.resize(container.buffer_data_uncompressed.size() + 65536);
		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);

		libzpaq::compress(&in, &out, &command[0]);

		const float fold = (float)container.buffer_data_uncompressed.size() / out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_data.data(),
                   container.buffer_data_uncompressed.data(),
                   container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();
			if(compress_strides){
				if(container.header.data_header.hasMixedStride())
					return(this->compressStrides(container, command));
				else return true;
			} else return true;
		}

		container.header.data_header.cLength            = out.buffer.size();
		container.header.data_header.controller.encoder = YON_ENCODE_ZPAQ;

		if(compress_strides){
			if(container.header.data_header.hasMixedStride())
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
			memcpy(container.buffer_data.data(),
                   container.buffer_data_uncompressed.data(),
                   container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.hasMixedStride())
					return(this->compressStrides(container));
				else return true;
			} else return true;
			//return true;
		}

		container.buffer_data.reset();
		container.buffer_data.resize(container.buffer_data_uncompressed.size() + 65536);
		ZpaqWrapperIn  in(container.buffer_data_uncompressed);
		ZpaqWrapperOut out(container.buffer_data);
		libzpaq::compress(&in, &out, "x0.3ci1");

		const float fold = (float)container.buffer_data_uncompressed.size() / out.buffer.size();
		if(fold < MIN_COMPRESSION_FOLD){
			memcpy(container.buffer_data.data(),
                   container.buffer_data_uncompressed.data(),
                   container.buffer_data_uncompressed.size());
			container.header.data_header.controller.encoder = YON_ENCODE_NONE;
			container.buffer_data.n_chars                   = container.buffer_data_uncompressed.size();
			container.header.data_header.cLength            = container.buffer_data_uncompressed.size();

			if(compress_strides){
				if(container.header.data_header.hasMixedStride())
					return(this->compressStrides(container));
				else return true;
			} else return true;
		}

		container.header.data_header.cLength            = out.buffer.size();
		container.header.data_header.controller.encoder = YON_ENCODE_ZPAQ;

		if(compress_strides){
			if(container.header.data_header.hasMixedStride())
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

		container.buffer_strides.reset();
		ZpaqWrapperIn  in(container.buffer_strides_uncompressed);
		ZpaqWrapperOut out(container.buffer_strides);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, &command[0]);

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
			container.buffer_strides.n_chars                  = container.buffer_strides_uncompressed.size();
			container.header.stride_header.cLength            = container.buffer_strides_uncompressed.size();
			return true;
		}

		container.buffer_strides.reset();
		ZpaqWrapperIn  in(container.buffer_strides_uncompressed);
		ZpaqWrapperOut out(container.buffer_strides);
		out.buffer.resize(in.buffer.size() + 65536);
		libzpaq::compress(&in, &out, "x0.3ci1");

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

	const bool decompress(container_type& container){
		if(container.header.data_header.controller.encoder != YON_ENCODE_ZPAQ){
			return true;
		}

		container.buffer_data_uncompressed.reset();
		container.buffer_data_uncompressed.resize(container.header.data_header.uLength + 65536);
		ZpaqWrapperIn in(container.buffer_data);
		ZpaqWrapperOut out(container.buffer_data_uncompressed);
		libzpaq::decompress(&in, &out);
		std::cerr << "zpaq decode: " << container.buffer_data.size() << "->" << container.buffer_data_uncompressed.size() << std::endl;

		assert(out.buffer.size() == container.header.data_header.uLength);
		assert(container.checkCRC(0));

		return true;
	}
	const bool decompressStrides(container_type& container){
		if(container.header.stride_header.controller.encoder != YON_ENCODE_ZPAQ){
			return true;
		}

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



#endif /* ALGORITHM_COMPRESSION_ZPAQ_CODEC_H_ */
