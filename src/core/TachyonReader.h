#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include "zstd.h"
#include "common/zstd_errors.h"
#include "BlockEntry.h"

namespace Tachyon{
namespace Core{

class TachyonReader{
	typedef TachyonReader self_type;
	typedef Core::BlockEntry block_entry_type;
	typedef IO::BasicBuffer buffer_type;

public:

	TachyonReader() : filesize(0){}
	TachyonReader(const std::string& filename) : input_file(filename), filesize(0){}
	~TachyonReader(){}

	bool open(void);
	bool open(const std::string& filename){
		if(filename.size() == 0){
			std::cerr << "no filename" << std::endl;
			return false;
		}
		this->stream.open(filename, std::ios::binary | std::ios::in | std::ios::ate);
		this->filesize = (U64)this->stream.tellg();
		if(!this->stream.good()){
			std::cerr << "failed to read file" << std::endl;
			return false;
		}
		this->stream.seekg(0);
		if(!this->stream.good()){
			std::cerr << "failed to rewrind" << std::endl;
			return false;
		}
		return true;
	}

	bool nextBlock(){
		if(!this->stream.good()){
			std::cerr << "faulty stream" << std::endl;
			return false;
		}

		if((U64)this->stream.tellg() == this->filesize){
			std::cerr << "eof all done" << std::endl;
			return false;
		}

		this->stream >> this->block;
		this->block.meta_hot_container.buffer_data_uncompressed.resize(this->block.meta_hot_container.header.uLength + 16536);
		if(this->block.meta_hot_container.header.controller.encoder == Core::ENCODE_ZSTD){
			int ret = ZSTD_decompress(this->block.meta_hot_container.buffer_data_uncompressed.data,
									  this->block.meta_hot_container.buffer_data_uncompressed.capacity(),
									  this->block.meta_hot_container.buffer_data.data,
									  this->block.meta_hot_container.buffer_data.pointer);

			assert(ret >= 0);
			this->block.meta_cold_container.buffer_data_uncompressed.pointer = ret;
			//std::cerr << "de: " << ret << " expected: " << this->block.meta_hot_container.header.uLength << std::endl;
			assert((U32)ret == this->block.meta_hot_container.header.uLength);
		}

		if(this->block.gt_rle_container.header.controller.encoder == Core::ENCODE_ZSTD){
			this->block.gt_rle_container.buffer_data_uncompressed.resize(this->block.gt_rle_container.header.uLength + 16536);
			int retRLE = ZSTD_decompress(this->block.gt_rle_container.buffer_data_uncompressed.data,
									  this->block.gt_rle_container.buffer_data_uncompressed.capacity(),
									  this->block.gt_rle_container.buffer_data.data,
									  this->block.gt_rle_container.buffer_data.pointer);
			assert(retRLE >= 0);
			this->block.gt_rle_container.buffer_data_uncompressed.pointer = retRLE;
			assert((U32)retRLE == this->block.gt_rle_container.header.uLength);
		}

		const Core::EntryHotMeta<U16>* const meta = reinterpret_cast<const Core::EntryHotMeta<U16>* const>(this->block.meta_hot_container.buffer_data_uncompressed.data);
		for(U32 i = 0; i < this->block.index_entry.n_variants; ++i)
			std::cout << meta[i] << '\n';

		return true;
	}
	bool seekBlock(const U32& b);

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	block_entry_type block;
	buffer_type buffer;
};

}
}

#endif /* CORE_TACHYONREADER_H_ */
