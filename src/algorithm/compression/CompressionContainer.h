#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../../io/compression/TGZFController.h"

namespace Tachyon{
namespace Compression{

inline bool bytePreprocessBits(const char* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data);
	BYTE* target[8];
	const U32 s = size/8;
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &dest[s*i];

	U32 k = 0; U32 p = 0;
	for(U32 i = 0; i < size; ++i, ++k){
		for(U32 j = 0; j < 8; ++j)
			target[j][p] |= ((d[i] & (1 << j)) >> j) << k;

		if(i % 7 == 0){ k = 0; ++p; }
	}

	return true;
}

class CompressionContainer{
private:
	typedef CompressionContainer self_type;
protected:
	typedef Core::StreamContainer stream_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::PermutationManager permutation_type;

public:
	CompressionContainer(){}
	virtual ~CompressionContainer(){}
	virtual const bool encode(permutation_type& manager) =0;
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
		stream.header.controller.encoder = Core::ENCODE_DEFLATE;
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
		stream.header.controller.encoder = Core::ENCODE_DEFLATE;
		std::cerr << "DEFLATE: " << stream.buffer_strides.pointer << '\t' << this->controller.buffer.size() << std::endl;
		memcpy(stream.buffer_strides.data, this->controller.buffer.data, this->controller.buffer.pointer);
		stream.buffer_strides.pointer = this->controller.buffer.pointer;
		this->controller.Clear();
		return true;
	}

	const bool encode(permutation_type& manager){
		this->buffer.reset();
		this->buffer.resize(manager.n_samples*sizeof(U32));
		// First
		const BYTE w = ceil(log2(manager.n_samples+1) / 8);
		if(w == 1){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (BYTE)manager[i];
		} else if(w == 2){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U16)manager[i];
			memset(manager.PPA.data, 0, this->buffer.pointer);
			manager.PPA.pointer = this->buffer.pointer;
			const U32 block_size = this->buffer.pointer / w;
			std::cerr << "repacking to: " << block_size << " @ " << (int)w << std::endl;
			bytePreprocessBits(&this->buffer.data[0], manager.PPA.pointer, &manager.PPA.data[0]);

			//for(U32 i = 0; i < manager.PPA.pointer; ++i)
			//	std::cerr << std::bitset<8>(manager.PPA.data[i]);
			//std::cerr << std::endl;

		} else if(w == 3 || w == 4){
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U32)manager[i];
		} else {
			for(U32 i = 0; i < manager.n_samples; ++i) this->buffer += (U64)manager[i];
		}
		//memcpy(manager.PPA.data, this->buffer.data, this->buffer.pointer);
		//manager.PPA.pointer = this->buffer.pointer;



		this->controller.Deflate(manager.PPA);
		manager.c_length = this->controller.buffer.size();
		//stream.header.controller.encoder = Core::ENCODE_DEFLATE;
		std::cerr << "DEFLATE PPA: " << manager.PPA.pointer << '\t' << this->controller.buffer.size() << std::endl;
		memcpy(manager.PPA.data, this->controller.buffer.data, this->controller.buffer.pointer);
		manager.PPA.pointer = this->controller.buffer.pointer;
		this->controller.Clear();
		return true;
	}

	const bool decode(stream_type& stream){ return false; }

private:
	controller_type controller;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
