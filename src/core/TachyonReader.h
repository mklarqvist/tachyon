#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include "zstd.h"
#include "common/zstd_errors.h"
#include "BlockEntry.h"
#include "../algorithm/compression/CompressionManager.h"
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
	typedef Compression::CompressionManager codec_manager_type;
	typedef BlockEntrySettings settings_type;

public:
	TachyonReader() : filesize(0){}
	TachyonReader(const std::string& filename) : input_file(filename), filesize(0){}
	~TachyonReader(){ }

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
		const U64 return_pos = this->stream.tellg();

		this->stream.seekg(this->filesize - 32);
		BYTE eof_data[32];
		Helpers::HexToBytes(Constants::TACHYON_FILE_EOF, &eof_data[0]);
		BYTE eof_match[32];
		this->stream.read((char*)&eof_match[0], 32);
		for(U32 i = 0; i < 32; ++i){
			if(eof_data[i] != eof_match[i]){
				std::cerr << "File is truncated!" << std::endl;
				return false;
			}
		}
		this->stream.seekg(this->filesize - 32 - sizeof(U64));
		this->stream.read((char*)reinterpret_cast<char*>(&this->l_data), sizeof(U64));
		std::cerr << this->l_data << std::endl;
		this->stream.seekg(return_pos);

		return true;
	}

	bool nextBlock(){
		if(!this->stream.good()){
			std::cerr << "faulty stream" << std::endl;
			return false;
		}

		if((U64)this->stream.tellg() == this->l_data){
			std::cerr << "eof all done" << std::endl;
			return false;
		}

		std::cerr << "before read" << std::endl;
		this->block.read(stream, this->settings);

		std::cerr << "after read" << std::endl;

		// Phase 1: Decode data
		// Todo:API: decode available data
		if(!this->codec_manager.decompress(this->block.meta_hot_container)){ std::cerr << "failed to decompress!" << std::endl; }
		if(!this->codec_manager.decompress(this->block.meta_cold_container)){ std::cerr << "failed to decompress!" << std::endl; }
		if(!this->codec_manager.decompress(this->block.gt_rle_container)){ std::cerr << "failed to decompress!" << std::endl; }
		if(!this->codec_manager.decompress(this->block.gt_simple_container)){ std::cerr << "failed to decompress!" << std::endl; }
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			//std::cerr << "checking: " << i << std::endl;
			if(!this->codec_manager.decompress(this->block.info_containers[i])){ std::cerr << "failed to decompress!" << std::endl; }
		}
		for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
			if(!this->codec_manager.decompress(this->block.format_containers[i])){ std::cerr << "failed to decompress!" << std::endl; }
		}

		return true;
	}

	bool seekBlock(const U32& b);

	bool toVCFPartial(void){
		// Phase 1 construct iterators
		Iterator::MetaIterator it;
		Iterator::ContainerIterator* info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];

		if(this->settings.loadMetaHot && !this->settings.loadMetaCold){
			it.set(this->block.meta_hot_container);
		} else if(this->settings.loadMetaHot && this->settings.loadMetaCold){
			it.set(this->block.meta_hot_container, this->block.meta_cold_container);
		} else if(this->settings.loadMetaCold && !this->settings.loadMetaHot){
			std::cerr << "illegal to load only cold meta" << std::endl;
			return false;
		}

		if(this->settings.load_info_ID_loaded.size() > 0 && this->settings.loadMetaHot == false){
			std::cerr << "not legal" << std::endl;
			return false;
		}

		// Setup containers
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i)
			info_iterators[i].setup(this->block.info_containers[i]);

		// Phase 2 perform iterations
		for(U32 i = 0; i < this->block.index_entry.n_variants; ++i){
			const Core::MetaEntry& m = it.current();
			std::cout.write(&this->header.getContig(this->block.index_entry.contigID).name[0], this->header.getContig(this->block.index_entry.contigID).name.size()) << '\t';
			std::cout << this->block.index_entry.minPosition + m.hot->position + 1 << '\t';

			// If we have cold meta
			if(this->settings.loadMetaCold){
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


				if(isnan(m.cold.QUAL)) std::cout << "\t.\t";
				else std::cout << '\t' << m.cold.QUAL << '\t';
			} else {
				std::cout << "\t.\t.\t.\t";
			}

			// Filter streams
			if(this->block.index_entry.n_filter_streams == 0){
				std::cout << ".\t";
			} else {
				for(U32 k = 0; k < this->block.index_entry.n_filter_streams; ++k){
					// Check if field is set
					if(this->block.index_entry.filter_bit_vectors[m.hot->FILTER_map_ID][k]){
					// Lookup what that field is
						std::cout.write(&this->header.getEntry(this->header.mapTable[this->block.index_entry.filter_offsets[k].key]).ID[0], this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID.size()) << '\t';
					}
				}
			}

			// Cycle over streams that are set in the given bit-vector
			U32 set = 0;
			const Index::IndexBlockEntryBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[m.hot->INFO_map_ID];
			//bool seen11 = false;
			for(U32 k = 0; k < this->settings.load_info_ID_loaded.size(); ++k){
				// If this key is set
				//std::cerr << k << '\t' << this->settings.load_info_ID_loaded[k].key << '\t' << this->settings.load_info_ID_loaded[k].target_stream << '\t' << this->settings.load_info_ID_loaded[k].target_stream_local << " validate: " << this->block.index_entry.info_offsets[this->settings.load_info_ID_loaded[k].target_stream_local].key << std::endl;

				if(target_info_vector[this->settings.load_info_ID_loaded[k].target_stream_local]){
					//std::cerr << "match: " << k << '\t' << this->settings.load_info_ID_loaded[k].target_stream_local << '\t' << this->settings.load_info_ID_loaded[k].offset->key << std::endl;
					if(set++ != 0)
						std::cout.put(';');

					info_iterators[this->settings.load_info_ID_loaded[k].iterator_index].toString(
							std::cout,
							this->header.getEntry(this->header.mapTable[this->settings.load_info_ID_loaded[k].offset->key]).ID);

					//if(set + 1 != this->settings.load_info_ID_loaded.size())


					++info_iterators[this->settings.load_info_ID_loaded[k].iterator_index];
					//if(this->settings.load_info_ID_loaded[k].key == 11) exit(1);
					//++set;
				}
			}
			//if(seen11) exit(1);
			std::cout.put('\n');

			++it;
		}
		std::cout.flush();

		delete [] info_iterators;
		return true;
	}

	bool toVCF(void){
		// Phase 1 construct iterators
		Iterator::MetaIterator it(this->block.meta_hot_container, this->block.meta_cold_container);
		Iterator::ContainerIterator* info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];

		// Setup containers
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			//std::cerr << this->header.entries[this->block.index_entry.info_offsets[i].key].ID << "\t" << this->block.info_containers[i].header.cLength << std::endl;
			info_iterators[i].setup(this->block.info_containers[i]);
		}

		// Todo: format

		// Phase 2 perform iterations
		for(U32 i = 0; i < this->block.index_entry.n_variants; ++i){
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

			if(isnan(m.cold.QUAL)) std::cout << "\t.\t";
			else std::cout << '\t' << m.cold.QUAL << '\t';

			if(this->block.index_entry.n_filter_streams == 0){
				std::cout << ".\t";
			} else {
				for(U32 k = 0; k < this->block.index_entry.n_filter_streams; ++k){
					// Check if field is set
					if(this->block.index_entry.filter_bit_vectors[m.hot->FILTER_map_ID][k]){
					// Lookup what that field is
						std::cout.write(&this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID[0], this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID.size()) << '\t';
					}
				}
			}

			// Cycle over streams that are set in the given bit-vector
			U32 set = 0;
			const Index::IndexBlockEntryBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[m.hot->INFO_map_ID];
			const U32 n_keys = target_info_vector.n_keys;
			const U32* const firstKey = &target_info_vector.keys[0];

			//for(U32 i = 0; i < this->header.n_entries; ++i){
			//	std::cerr << i << '\t' << this->header.entries[i].ID << std::endl;
			//}

			//std::cerr << "\nSet membership: " << m.hot->INFO_map_ID << std::endl;
			//for(U32 k = 0; k < n_keys; ++k){
			//	std::cerr << firstKey[k] << "->" << this->block.index_entry.info_offsets[firstKey[k]].key << '\t' << this->header.entries[this->block.index_entry.info_offsets[firstKey[k]].key].ID << std::endl;
			//}
			//std::cerr << std::endl;

			for(U32 k = 0; k < n_keys; ++k){
				// Check if field is set
				const U32& current_key = firstKey[k];
				if(target_info_vector[current_key]){
					info_iterators[current_key].toString(std::cout, this->header.entries[this->header.mapTable[this->block.index_entry.info_offsets[firstKey[k]].key]].ID);

					if(set + 1 != target_info_vector.n_keys)
						std::cout.put(';');

					++info_iterators[current_key];
					++set;
				}
			}
			std::cout << '\n';
			++it;
		}
		std::cout.flush();

		delete [] info_iterators;
		return true;
	}

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	U64 l_data;

	settings_type settings;
	block_entry_type block;
	header_type header;

	codec_manager_type codec_manager;
};

}
}

#endif /* CORE_TACHYONREADER_H_ */
