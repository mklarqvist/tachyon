#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"

namespace Tachyon{
namespace Compression{

class CompressionContainer{
private:
	typedef CompressionContainer self_type;
protected:
	typedef Core::StreamContainer stream_type;
	typedef IO::BasicBuffer buffer_type;

public:
	CompressionContainer(){}
	virtual ~CompressionContainer(){}
	virtual const bool encode(stream_type& stream) =0;
	virtual const bool encodeStrides(stream_type& stream) =0;
	virtual const bool decode(stream_type& stream) =0;

protected:
	buffer_type buffer;
};

// GZIP of delfate codec
class DeflateCodec : public CompressionContainer{
private:
	typedef DeflateCodec self_type;
	typedef IO::TGZFController controller_type;

public:
	DeflateCodec(){}
	~DeflateCodec(){ }

	const bool encode(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_data.pointer);
		this->controller.Deflate(stream.buffer_data);
		stream.header.uLength = stream.buffer_data.pointer;
		stream.header.cLength = this->controller.buffer.size();
		std::cerr << "DEFLATE: " << stream.buffer_data.pointer << '\t' << this->controller.buffer.size() << std::endl;
		memcpy(stream.buffer_data.data, this->controller.buffer.data, this->controller.buffer.pointer);
		stream.buffer_data.pointer = this->controller.buffer.pointer;
		this->controller.Clear();
		return true;
	}

	const bool encodeStrides(stream_type& stream){
		this->buffer.reset();
		this->buffer.resize(stream.buffer_strides.pointer);
		this->controller.Deflate(stream.buffer_strides);
		stream.header.uLength = stream.buffer_strides.pointer;
		stream.header.cLength = this->controller.buffer.size();
		std::cerr << "DEFLATE: " << stream.buffer_strides.pointer << '\t' << this->controller.buffer.size() << std::endl;
		memcpy(stream.buffer_strides.data, this->controller.buffer.data, this->controller.buffer.pointer);
		stream.buffer_strides.pointer = this->controller.buffer.pointer;
		this->controller.Clear();
		return true;
	}

	bool encrypt(stream_type& stream){
		const std::string key_hex = "FC6F9E22E83EFF8C8120CF233E05CB2EFCC20EBB1212DF44AF919184595A5355";
		const std::string iv_hex = "E3C2AA42EE6F0DBB7556D1D38290DE7A";

		// Todo: fix
		uint8_t key[32]; memset(key, 0x00, 32);
		uint8_t iv[16];  memset(iv, 0x00, 16);
		Helpers::HexToBytes(key_hex, &key[0]);
		Helpers::HexToBytes(iv_hex, &iv[0]);
		//AES_CBC_encrypt_buffer(reinterpret_cast<uint8_t*>(this->buffer.data), reinterpret_cast<uint8_t*>(stream.buffer_data.data), stream.buffer_data.pointer, &key[0], &iv[0]);
		//AES_ECB_encrypt(reinterpret_cast<uint8_t*>(stream.buffer_data.data),&key[0],reinterpret_cast<uint8_t*>(this->buffer.data),stream.buffer_data.pointer);
		return(true);
	}

	const bool decode(stream_type& stream){ return false; }

private:
	controller_type controller;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
