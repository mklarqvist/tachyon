#ifndef ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_
#define ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_

#include "libzpaq.h"
#include "io/basic_buffer.h"

namespace tachyon{
namespace algorithm{

class ZpaqWrapperIn: public libzpaq::Reader {
private:
	typedef io::BasicBuffer buffer_type;

public:
	ZpaqWrapperIn(buffer_type& buffer) : iterator_pos(0), buffer(buffer){}
	~ZpaqWrapperIn(){ }
	inline int get(){
		if(this->iterator_pos + 1 == this->buffer.size()) return(-1); // eof
		return((BYTE)this->buffer[this->iterator_pos++]);
		getchar();
	}  // returns byte 0..255 or -1 at EOF

	inline void reset(void){
		this->buffer.reset();
		this->iterator_pos = 0;
	}

public:
	size_t       iterator_pos;
	buffer_type& buffer;
 };

 class ZpaqWrapperOut: public libzpaq::Writer {
 private:
 	typedef io::BasicBuffer buffer_type;

 public:
	ZpaqWrapperOut(buffer_type& buffer) : buffer(buffer){}
	~ZpaqWrapperOut(){ }
	inline void put(int c){ this->buffer += (BYTE)c; }  // writes 1 byte 0..255
	//inline void write(const char* buf, int n){ this->buffer.Add(buf, n); }
	inline void reset(void){ this->buffer.reset(); }

 public:
	buffer_type& buffer;
};

}
}

#endif /* ALGORITHM_COMPRESSION_ZPAQ_WRAPPER_H_ */
