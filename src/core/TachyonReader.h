#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"
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

		// Reset
		this->block.clear();

		if(!this->block.read(stream, this->settings))
			return false;

		// Phase 1: Decode data
		if(!this->codec_manager.decompress(this->block))
			return false;

		return true;
	}

	bool seekBlock(const U32& b);

	bool toVCFString(std::ostream& stream = std::cout){
		// Phase 1 construct iterators
		Iterator::MetaIterator* it = this->block.getMetaIterator(); // factory

		// Setup containers
		// INFO
		Iterator::ContainerIterator* info_iterators;
		if(this->settings.loadInfoAll){
			info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];

			for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
				info_iterators[i].setup(this->block.info_containers[i]);
			}
		}

		// FORMAT
		Iterator::ContainerIterator* format_iterators;
		if(this->settings.loadFormatAll){
			format_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_format_streams];
			for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
				format_iterators[i].setup(this->block.format_containers[i]);
			}
		}

		// Phase 2 perform iterations

		for(U32 i = 0; i < this->block.index_entry.size(); ++i){
			const Core::MetaEntry& m = (*it)[i];

			// Discussion: it is more attractive to have the
			// struct have knowledge of itself
			// it is not pretty to have to pass along
			// references to the header and index
			// every time
			m.toVCFString(stream,
					this->header,
					this->block.index_entry.contigID,
					this->block.index_entry.minPosition);

			if(this->block.index_entry.n_filter_streams == 0){
				stream << ".\t";
			} else {
				for(U32 k = 0; k < this->block.index_entry.n_filter_streams; ++k){
					// Check if field is set
					if(this->block.index_entry.filter_bit_vectors[m.filter_pattern_id][k]){
					// Lookup what that field is
						stream.write(&this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID[0], this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID.size()) << '\t';
					}
				}
			}

			U32 set = 0;
			if(this->settings.loadInfoAll){
				// Cycle over streams that are set in the given bit-vector

				const Index::IndexBlockEntryBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[m.info_pattern_id];
				const U32 n_keys = target_info_vector.n_keys;
				const U32* const firstKey = &target_info_vector.keys[0];

				for(U32 k = 0; k < n_keys; ++k){
					// Check if field is set
					const U32& current_key = firstKey[k];
					if(target_info_vector[current_key]){
						info_iterators[current_key].toString(stream, this->header.entries[this->header.mapTable[this->block.index_entry.info_offsets[firstKey[k]].key]].ID);

						if(set + 1 != target_info_vector.n_keys)
							stream.put(';');

						++info_iterators[current_key];
						++set;
					}
				}
			}

			if(this->settings.loadFormatAll){

				// Start FORMAT description field
				const Index::IndexBlockEntryBitvector& target_format_vector = this->block.index_entry.format_bit_vectors[m.format_pattern_id];
				const U32 n_keys_format = target_format_vector.n_keys;
				const U32* const firstKey_format = &target_format_vector.keys[0];

				stream.put('\t');
				if(n_keys_format){
					stream << this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[firstKey_format[0]].key]].ID;
					for(U32 j = 1; j < n_keys_format; ++j){
						stream.put(':');
						stream << this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[firstKey_format[j]].key]].ID;
					}
				} else {
					goto end_format;
				}
				stream.put('\t');

				// Format here
				// Has to admix iterators
				// Cycle over streams that are set in the given bit-vector
				set = 0;

				for(U32 s = 0; s < this->header.n_samples; ++s){
					for(U32 k = 0; k < n_keys_format; ++k){
						// Check if field is set
						const U32& current_key = firstKey_format[k];
						if(target_format_vector[current_key]){
							//std::cerr << '\t' << current_key << "->" << this->block.index_entry.format_offsets[current_key].key << ";" << this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[current_key].key]].ID << '\t' << format_iterators[current_key].data_iterator->n_entries << '\t' << (U32)format_iterators[current_key].get<BYTE>() << std::endl;
							if(this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[firstKey_format[k]].key]].ID == "GT"){
								//std::cerr << "isGT" << std::endl;
								//++format_iterators[current_key];
								++set;
								continue;
							}

							format_iterators[current_key].toString(stream);

							//std::cout << '@';

							if(set + 1 != target_format_vector.n_keys)
								stream.put(':');

							format_iterators[current_key].increment(false);
							++set;
						}
					}
					set = 0;
					stream.put('\t');
				}

				for(U32 k = 0; k < n_keys_format; ++k){
					if(this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[firstKey_format[k]].key]].ID == "GT"){
						//std::cerr << "isGT" << std::endl;
						//++format_iterators[current_key];

						continue;
					}
					// Check if field is set
					const U32& current_key = firstKey_format[k];
					format_iterators[current_key].incrementStride();
				}
			}

			end_format:
			stream << '\n';

		}
		stream.flush();

		delete it;
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
