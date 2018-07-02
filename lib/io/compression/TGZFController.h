#ifndef TGZFCONTROLLER_H_
#define TGZFCONTROLLER_H_

#include <fstream>

#include "../../support/helpers.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"
#include "../basic_buffer.h"
#include "GZFHeader.h"

namespace tachyon{
namespace io{

class TGZFController{
private:
	typedef TGZFController self_type;

protected:
	typedef io::BasicBuffer buffer_type;
	typedef TGZFHeader header_type;

public:
	TGZFController();
	TGZFController(const char* data, const U32 length);
	TGZFController(const U32 largest_block_size);
	~TGZFController();

	void Clear();
	bool Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;
	bool Inflate(buffer_type& input, buffer_type& output) const;
	bool InflateBlock(std::ifstream& stream, buffer_type& input);

	bool Deflate(const buffer_type& buffer);
	bool Deflate(buffer_type& meta, buffer_type& rle);
	bool Deflate(buffer_type& meta, buffer_type& meta_complex, buffer_type& rle);
	inline void setWindowSize(const S32& window){ this->bit_window = window; }
	inline void setCompression(const S32& compression){ this->compression_level = compression; }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(entry.buffer.buffer, entry.buffer.n_chars);
		return stream;
	}

private:
	bool __Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;

public:
	S32 compression_level;
	S32 bit_window;
	buffer_type buffer;
};

}
}

#endif /* TGZFCONTROLLER_H_ */
