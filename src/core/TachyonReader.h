#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include "zstd.h"
#include "common/zstd_errors.h"
#include "BlockEntry.h"
#include "../algorithm/compression/CompressionContainer.h"
#include "decorator/MetaHotDecorator.h"
#include "base/MetaCold.h"
#include "base/MetaEntry.h"

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

		if(this->block.meta_cold_container.header.controller.encoder == Core::ENCODE_ZSTD){
			this->zstd.decode(this->block.meta_cold_container);
		}

		if(this->block.meta_hot_container.header.controller.encoder == Core::ENCODE_ZSTD){
			this->zstd.decode(this->block.meta_hot_container);

			// Todo: MetaIterator
			Decorator::MetaHotDecorator d(this->block.meta_hot_container, this->block.index_entry.minPosition);
			//const Core::EntryHotMeta* hot = nullptr;
			//std::cerr << this->block.index_entry.maxPosition - this->block.index_entry.minPosition << " bp" << std::endl;
			//U32 prevpos = d[0].position;
			for(U32 i = 0; i < d.size(); ++i){
				//std::cout << d[i] << '\t' << this->block.index_entry.minPosition + d[i].position << '\t' << d[i].position - prevpos  << '\n';
				//const Core::MetaCold& cold = *reinterpret_cast<const Core::MetaCold* const>(&this->block.meta_cold_container.buffer_data_uncompressed[d[i].virtual_offset_cold_meta]);
				Core::MetaEntry m(&d[i], this->block.meta_cold_container);
				//std::cout << *m.hot << '\n';
				std::cout << this->block.index_entry.contigID << '\t';
				std::cout << this->block.index_entry.minPosition + m.hot->position + 1 << '\t';
				if(m.cold.n_ID == 0) std::cout.put('.');
				else std::cout.write(m.cold.ID, m.cold.n_ID);
				std::cout << '\t';
				if(m.hot->controller.biallelic && m.hot->controller.simple){
					std::cout << m.hot->ref_alt.getRef() << '\t' << m.hot->ref_alt.getAlt();
				}
				else {
					std::cout.write(m.cold.alleles[0].allele, m.cold.alleles[0].l_allele);
					std::cout << '\t';
					U32 i = 1;
					for(; i < m.cold.n_allele - 1; ++i){
						std::cout.write(m.cold.alleles[i].allele, m.cold.alleles[i].l_allele);
						std::cout.put(',');
					}
					std::cout.write(m.cold.alleles[i].allele, m.cold.alleles[i].l_allele);
				}
				std::cout << '\t' << m.cold.QUAL << '\t' << "PASS" << '\n';
			}
		}

		/*
		if(this->block.gt_rle_container.header.controller.encoder == Core::ENCODE_ZSTD){
			this->zstd.decode(this->block.gt_rle_container);
		}
		if(this->block.gt_simple_container.header.controller.encoder == Core::ENCODE_ZSTD){
			this->zstd.decode(this->block.gt_simple_container);
		}

		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			if(this->block.info_containers[i].header.controller.encoder == Core::ENCODE_ZSTD){
				this->zstd.decode(this->block.info_containers[i]);
			} else if(this->block.info_containers[i].header.controller.encoder == Core::ENCODE_NONE){
				std::cerr << "ENCODE_NONE | CRC check " << (this->block.info_containers[i].checkCRC(3) ? "PASS" : "FAIL") << std::endl;
			}
		}

		for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
			if(this->block.format_containers[i].header.controller.encoder == Core::ENCODE_ZSTD){
				this->zstd.decode(this->block.format_containers[i]);
			} else if(this->block.format_containers[i].header.controller.encoder == Core::ENCODE_NONE){
				std::cerr << "ENCODE_NONE | CRC check " << (this->block.format_containers[i].checkCRC(3) ? "PASS" : "FAIL") << std::endl;
			}
		}
		*/

		/*
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
		*/

		//const Core::EntryHotMeta<U16>* const meta = reinterpret_cast<const Core::EntryHotMeta<U16>* const>(this->block.meta_hot_container.buffer_data_uncompressed.data);
		//for(U32 i = 0; i < this->block.index_entry.n_variants; ++i)
		//	std::cout << meta[i] << '\n';

		return true;
	}
	bool seekBlock(const U32& b);

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	block_entry_type block;
	buffer_type buffer;
	Compression::ZSTDCodec zstd;
};

}
}

#endif /* CORE_TACHYONREADER_H_ */
