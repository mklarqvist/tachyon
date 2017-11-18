#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

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
	~DeflateCodec(){}

	const bool encode(stream_type& stream){
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

	const bool decode(stream_type& stream){ return false; }

private:
	controller_type controller;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
