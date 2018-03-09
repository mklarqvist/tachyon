#ifndef ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_
#define ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_

#include "libzpaq.h"
#include "../../io/basic_buffer.h"

namespace tachyon{
namespace algorithm{

class ZpaqWrapperIn: public libzpaq::Reader {
public:
	ZpaqWrapperIn(io::BasicBuffer& buffer) : iterator_pos(0), buffer(buffer){}
	~ZpaqWrapperIn(){ }
	inline int get() {
		//std::cerr << this->iterator_pos << "/" << this->buffer.size() << ":" << (int)this->buffer[this->iterator_pos] << ' ';
		if(this->iterator_pos + 1 == this->buffer.size()) return(-1); // eof
		return(*reinterpret_cast<const BYTE* const>(&this->buffer[this->iterator_pos++]));
	}  // returns byte 0..255 or -1 at EOF

	void reset(void){
		this->buffer.reset();
		this->iterator_pos = 0;
	}

	size_t           iterator_pos;
	io::BasicBuffer& buffer;
 };

 class ZpaqWrapperOut: public libzpaq::Writer {
 public:
	ZpaqWrapperOut(io::BasicBuffer& buffer) : buffer(buffer){}
	~ZpaqWrapperOut(){ }
	inline void put(int c) { this->buffer += (BYTE)c; }  // writes 1 byte 0..255
	void write(const char* buf, int n){  this->buffer.Add(buf, n); }
	void reset(void){ this->buffer.reset(); }

	io::BasicBuffer& buffer;
};

}
}

#endif /* ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_ */
