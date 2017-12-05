#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include "zstd.h"
#include "common/zstd_errors.h"
#include "BlockEntry.h"
#include "../algorithm/compression/CompressionContainer.h"
#include "iterator/MetaIterator.h"
#include "base/header/Header.h"
#include "iterator/ContainerIterator.h"

namespace Tachyon{
namespace Core{

class TachyonReader{
	typedef TachyonReader self_type;
	typedef Core::BlockEntry block_entry_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Header header_type;

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

		this->stream << this->header;
		if(!this->stream.good()){
			std::cerr << "failed to get header" << std::endl;
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
		}


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
				//std::cerr << "ENCODE_NONE | CRC check " << (this->block.info_containers[i].checkCRC(3) ? "PASS" : "FAIL") << std::endl;
				this->block.info_containers[i].buffer_data_uncompressed.resize(this->block.info_containers->buffer_data.pointer + 16536);
				memcpy(this->block.info_containers[i].buffer_data_uncompressed.data, this->block.info_containers[i].buffer_data.data, this->block.info_containers[i].buffer_data.pointer);
				this->block.info_containers[i].buffer_data_uncompressed.pointer = this->block.info_containers[i].buffer_data.pointer;
			}

			if(this->block.info_containers[i].header.controller.mixedStride){
				if(this->block.info_containers[i].header_stride.controller.encoder == Core::ENCODE_ZSTD){
					this->zstd.decodeStrides(this->block.info_containers[i]);
				} else if (this->block.info_containers[i].header_stride.controller.encoder == Core::ENCODE_NONE){
					this->block.info_containers[i].buffer_strides_uncompressed.resize(this->block.info_containers->buffer_strides.pointer + 16536);
					memcpy(this->block.info_containers[i].buffer_strides_uncompressed.data, this->block.info_containers[i].buffer_strides.data, this->block.info_containers[i].buffer_strides.pointer);
					this->block.info_containers[i].buffer_strides_uncompressed.pointer = this->block.info_containers[i].buffer_strides.pointer;
				}
			}

			/*
			std::cerr << Helpers::timestamp("LOG","ITERATOR") <<
							(this->block.info_containers[i].header.controller.uniform ? "UNIFORM" : "NON-UNIFORM") << '\t' <<
							this->block.info_containers[i].buffer_data_uncompressed.pointer << '\t' <<
							this->block.info_containers[i].header.controller.type << '\t' <<
							this->block.info_containers[i].header.controller.signedness << "\tcompressed: " <<
							this->block.info_containers[i].buffer_data.pointer << '\t' <<
							this->header.getEntry(this->block.index_entry.info_offsets[i].key).ID << std::endl;
			*/

			//std::cerr << (this->block.info_containers[i].header.controller.type == Core::TYPE_BOOLEAN) << std::endl;
			//assert(this->block.info_containers[i].buffer_data_uncompressed.pointer > 0);
		}


		for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
			if(this->block.format_containers[i].header.controller.encoder == Core::ENCODE_ZSTD){
				this->zstd.decode(this->block.format_containers[i]);
			} else if(this->block.format_containers[i].header.controller.encoder == Core::ENCODE_NONE){
				//std::cerr << "ENCODE_NONE | CRC check " << (this->block.format_containers[i].checkCRC(3) ? "PASS" : "FAIL") << std::endl;
				this->block.format_containers[i].buffer_data_uncompressed.resize(this->block.format_containers->buffer_data.pointer + 16536);
				memcpy(this->block.format_containers[i].buffer_data_uncompressed.data, this->block.format_containers[i].buffer_data.data, this->block.format_containers[i].buffer_data.pointer);
				this->block.format_containers[i].buffer_data_uncompressed.pointer = this->block.format_containers[i].buffer_data.pointer;

			}
		}

		//std::cout << "first\n\n" << std::endl;

		Iterator::MetaIterator it(this->block.meta_hot_container, this->block.meta_cold_container);
		Iterator::ContainerIterator* info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			info_iterators[i].setup(this->block.info_containers[i]);
			std::cerr << i << "->" << this->header.getEntry(this->block.index_entry.info_offsets[i].key).ID << '\t' <<
					info_iterators[i].data_iterator->n_entries << '\t';
			if(info_iterators[i].stride_iterator == nullptr){
				std::cerr << "fixed : " << info_iterators[i].container->header.stride << std::endl;
			} else std::cerr << "variable: " << info_iterators[i].stride_iterator->n_entries << std::endl;
		}

		//std::cerr << it.size() << std::endl;
		for(U32 i = 0; i < it.size(); ++i){
			const Core::MetaEntry& m = it.current();
			std::cout.write(&this->header.getContig(this->block.index_entry.contigID).name[0], this->header.getContig(this->block.index_entry.contigID).name.size()) << '\t';
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
				U16 j = 1;
				for(; j < m.cold.n_allele - 1; ++j){
					std::cout.write(m.cold.alleles[j].allele, m.cold.alleles[j].l_allele);
					std::cout.put(',');
				}
				std::cout.write(m.cold.alleles[j].allele, m.cold.alleles[j].l_allele);
			}

			std::cout << '\t' << m.cold.QUAL << '\t';
			for(U32 k = 0; k < this->block.index_entry.n_filter_streams; ++k){
				// Check if field is set
				if(this->block.index_entry.filter_bit_vectors[m.hot->FILTER_map_ID][k]){
				// Lookup what that field is
					std::cout.write(&this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID[0], this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID.size()) << '\t';
				}
			}

			U32 set = 0;
			const Index::IndexBlockEntryBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[m.hot->INFO_map_ID];
			for(U32 k = 0; k < this->block.index_entry.n_info_streams; ++k){
				// Todo: we need to know the pattern order
				//if(k == 11) continue;
				//if(k == 16) continue;

				// Check if field is set
				if(target_info_vector[k]){
					//std::cout << this->block.index_entry.info_bit_vectors[m.hot->INFO_map_ID][k] << " is set: " << this->header.entries[k].ID << std::endl;
					// Lookup what that field is
					/*
					if( k > 4 ){
						//std::cout << this->header.getEntry(this->block.index_entry.info_offsets[k].key).ID << "=" << 0;
						std::cerr << "TYPE: " << (int)info_iterators[k].container->header.controller.type << std::endl;
						std::cerr << "OUTPUT: " << info_iterators[k].toString(std::cout, this->header.getEntry(this->block.index_entry.info_offsets[k].key).ID) << std::endl;
						//exit(1);
					} else {

					info_iterators[k].toString(std::cout, this->header.getEntry(this->block.index_entry.info_offsets[k].key).ID);
					}
					*/
					//if(k != 18)
						info_iterators[k].toString(std::cout, this->header.getEntry(this->block.index_entry.info_offsets[k].key).ID);
					/*
						else {
						std::cerr << k << " TYPE: " << (int)info_iterators[k].container->header.controller.type << std::endl;
						info_iterators[k].toString(std::cerr, this->header.getEntry(this->block.index_entry.info_offsets[k].key).ID);
						std::cerr << std::endl;
						std::cerr << info_iterators[k].container->buffer_data_uncompressed.data << std::endl;
						std::cerr << info_iterators[k].container->buffer_data.data << std::endl;
					}
					*/

					if(set + 1 != target_info_vector.fields_set)
						std::cout.put(';');

					++info_iterators[k];
					++set;
				}


			}
			std::cout << '\n';
			++it;
		}
		//exit(1);

		/*
		Iterator::ContainerIterator cit;
		cit(this->block.info_containers[0]);
		std::cerr << "type for iterator: " << cit.data_iterator->header.controller.type << std::endl;
		std::cerr << "mixed? : " << cit.data_iterator->header.controller.mixedStride << std::endl;
		std::cerr << "uniform? : " << cit.data_iterator->header.controller.uniform << std::endl;

		Iterator::ContainerIteratorData<U16> re_cit(cit.container->buffer_data_uncompressed, cit.container->header);

		if(this->block.info_containers[0].header.controller.mixedStride){
			std::cerr << "type for stride: " << cit.stride_iterator->header.controller.type << std::endl;

			Iterator::ContainerIteratorStride<BYTE> s_cit(cit.container->buffer_strides_uncompressed, cit.container->header_stride);
			std::cerr << this->header.entries[this->block.index_entry.info_offsets->key].ID << '\t' << re_cit.buffer.pointer << std::endl;

			for(U32 i = 0; i < s_cit.n_entries; ++i){
				for(U32 j = 0; j < s_cit.current(); ++j){
					std::cerr << re_cit.current() << ',';
					++re_cit;
				}
				std::cerr << std::endl;
				++s_cit;
			}
		}
		*/
		//std::cout << "last\n" << std::endl;

		//if((this->block.index_entry.minPosition + it.first().hot->position + 1 ) >= 118222)
		//	exit(1);

		delete [] info_iterators;
		return true;
	}
	bool seekBlock(const U32& b);

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	block_entry_type block;
	header_type header;
	buffer_type buffer;
	Compression::ZSTDCodec zstd;
};

}
}

#endif /* CORE_TACHYONREADER_H_ */
