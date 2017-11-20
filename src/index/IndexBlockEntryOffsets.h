#ifndef INDEX_INDEXBLOCKENTRYOFFSETS_H_
#define INDEX_INDEXBLOCKENTRYOFFSETS_H_

#include "../core/base/StreamContainerHeaderController.h"
#include "../core/base/StreamContainerHeader.h"

namespace Tachyon{
namespace Index{

struct IndexBlockEntryOffsets{
	typedef IndexBlockEntryOffsets self_type;
	typedef Core::StreamContainerHeader header_type;
	typedef Core::StreamContainerHeaderStride header_stride_type;

public:
	IndexBlockEntryOffsets(void) : key(0){}
	IndexBlockEntryOffsets(const U32& key, const header_type& h) : key(key), header(h){}
	IndexBlockEntryOffsets(const U32& key, const header_type& h, const header_stride_type& s) : key(key), header(h), header_stride(s){}
	~IndexBlockEntryOffsets(void){}

	bool update(const U32& key, const header_type& h){
		this->key = key;
		this->header = h;
		return true;
	}

	bool update(const U32& key, const header_type& h, const header_stride_type& s){
		this->key = key;
		this->header = h;
		this->header_stride = s;
		return true;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.key), sizeof(U32));
		stream << entry.header;
		if(entry.header.controller.mixedStride)
			stream << entry.header_stride;

		return(stream);
	}

public:
	U32 key;
	header_type header;
	header_stride_type header_stride;
};

}
}

#endif /* INDEX_INDEXBLOCKENTRYOFFSETS_H_ */
