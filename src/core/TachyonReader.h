#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"
#include "Block.h"
#include "../algorithm/compression/CompressionManager.h"
#include "iterator/MetaIterator.h"
#include "base/header/Header.h"
#include "iterator/ContainerIterator.h"
#include "../index/SortedIndex.h"
#include "iterator/GenotypeIterator.h"
#include "../algorithm/Timer.h"

namespace Tachyon{

class TachyonReader{
	typedef TachyonReader self_type;
	typedef Core::Block block_entry_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::Header header_type;
	typedef Compression::CompressionManager codec_manager_type;
	typedef Core::BlockEntrySettings settings_type;
	typedef Index::SortedIndex index_type;

public:
	TachyonReader() : filesize(0), n_internal_buffers(0), internal_buffers(nullptr){}
	TachyonReader(const std::string& filename) : input_file(filename), filesize(0), n_internal_buffers(0), internal_buffers(nullptr){}
	~TachyonReader(){
		for(U32 i = 0; i < this->n_internal_buffers; ++i)
			this->internal_buffers[i].deleteAll();

		delete [] this->internal_buffers;

	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void){
		if(this->input_file.size() == 0){
			std::cerr << "no filename" << std::endl;
			return false;
		}
		this->stream.open(this->input_file, std::ios::binary | std::ios::in | std::ios::ate);
		this->filesize = (U64)this->stream.tellg();
		if(!this->stream.good()){
			std::cerr << "failed to read file" << std::endl;
			return false;
		}
		this->stream.seekg(0);
		if(!this->stream.good()){
			std::cerr << "failed to rewind" << std::endl;
			return false;
		}

		// Load header
		this->stream << this->header;
		if(!this->stream.good()){
			std::cerr << "failed to get header" << std::endl;
			return false;
		}
		// Keep track of start position
		const U64 return_pos = this->stream.tellg();

		// Seek to EOF and make check
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
		// Seek back to find start of index
		// Seek to that start of index
		// Load index
		// Seek back to start of the file
		this->stream.seekg(this->filesize - 32 - sizeof(U64));
		this->stream.read((char*)reinterpret_cast<char*>(&this->l_data), sizeof(U64));
		this->stream.seekg(this->l_data);
		this->stream >> this->index;
		this->stream.seekg(return_pos);

		return(this->stream.good());
	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(const std::string& filename){
		this->input_file = filename;
		return(this->open());
	}

	/**<
	 * Overloaded operator for blocks. Useful when
	 * looping over a range of blocks. This happens
	 * frequently in parallel programs.
	 * @param index Block index value in range [0..n_blocks)
	 * @return      Returns TRUE if operation was successful or FALSE otherwise
	 */
	bool operator[](const U32 index);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool getNextBlock(){
		// If the stream is faulty then return
		if(!this->stream.good()){
			std::cerr << "faulty stream" << std::endl;
			return false;
		}

		// If the current position is the EOF then
		// exit the function
		if((U64)this->stream.tellg() == this->l_data)
			return false;

		// Reset and re-use
		this->block.clear();

		// Attempts to read a YON block with the provided
		// settings
		if(!this->block.read(stream, this->settings))
			return false;

		// Internally decompress available data
		if(!this->codec_manager.decompress(this->block))
			return false;

		// All passed
		return true;
	}

	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param b
	 * @return
	 */
	bool seekBlock(const U32& b);

	/**<
	 * Move all internal pointers up to the next available
	 * variant. If the next variant is in another block then
	 * load that block by invoking nextBlock().
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool nextVariant(void);

	bool toVCFStringFast(std::ostream& stream = std::cout){
		Iterator::MetaIterator* it = this->block.getMetaIterator(); // factory

		// Setup containers
		// INFO
		Iterator::ContainerIterator* info_iterators = nullptr;
		if(this->settings.loadInfoAll){
			info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];

			for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
				info_iterators[i].setup(this->block.info_containers[i]);
			}
		}

		// FORMAT
		Iterator::ContainerIterator* format_iterators = nullptr;
		if(this->settings.loadFormatAll){
			format_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_format_streams];
			for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
				format_iterators[i].setup(this->block.format_containers[i]);
			}
		}

		// One buffer per variant
		// Faster output?
		//std::cerr << this->block.size() << '/' << this->n_internal_buffers << std::endl;
		if(this->block.size() > this->n_internal_buffers){
			//std::cerr << "here: " << this->n_internal_buffers << std::endl;
			for(U32 i = 0; i < this->n_internal_buffers; ++i)
				this->internal_buffers[i].deleteAll();
			delete [] this->internal_buffers;

			this->n_internal_buffers = this->block.size() * 5;
			this->internal_buffers = new buffer_type[5*this->block.size()];
			for(U32 i = 0; i < 5*this->block.size(); ++i){
				this->internal_buffers[i].resize(5012);
				//std::cerr << i << '\t' << this->internal_buffers[i].capacity() << std::endl;
			}
		}

		// Base
		for(U32 i = 0; i < this->block.size(); ++i){
			const Core::MetaEntry& m = (*it)[i];

			m.toVCFString(this->internal_buffers[i],
                          this->header,
                          this->block.index_entry.contigID,
                          this->block.index_entry.minPosition);
		}

		// Memory layout:
		// Cycle over patterns
		// Cycle over entries
		// Cycle over each key in pattern
		// If key available: add to buffer
		if(this->settings.loadInfoAll){
			const U32 n_patterns = this->block.index_entry.n_info_patterns;
			for(U32 p = 0; p < n_patterns; ++p){
				// Cycle over streams that are set in the given bit-vector
				const Index::BlockIndexBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[p];
				const U32 n_keys = target_info_vector.n_keys;
				const U32* const keys = &target_info_vector.keys[0];

				for(U32 k = 0; k < n_keys; ++k){
					const U32& current_key = keys[k];

					for(U32 i = 0; i < this->block.size(); ++i){
						const Core::MetaEntry& m = (*it)[i];
						if(m.getInfoPatternID() != p)
							continue;

						info_iterators[current_key].toString(this->internal_buffers[i], this->header.entries[this->header.mapTable[this->block.index_entry.info_offsets[current_key].key]].ID);
						if(k + 1 != n_keys)
							this->internal_buffers[i] += ';';

						++info_iterators[current_key];
					}
				}
			}
		}

		for(U32 i = 0; i < this->block.size(); ++i){
			stream.write(this->internal_buffers[i].data, this->internal_buffers[i].pointer);
			stream.put('\n');
			this->internal_buffers[i].reset();
		}

		delete it;
		delete [] info_iterators;
		delete [] format_iterators;
		return true;
	}

	bool toVCFString(std::ostream& stream = std::cout){
		//Core::HeaderMapEntry* entry = nullptr;
		/*
		if(this->header.getEntry("GT", entry)){
			std::cerr << "GT@" << entry->ID << '\t' << entry->IDX << '\t' << (int)entry->TYPE << std::endl;
		}
		if(this->header.getEntry("PASS", entry)){
			std::cerr << "PASS@" << entry->ID << '\t' << entry->IDX << '\t' << (int)entry->TYPE << std::endl;
		}
		if(this->header.getEntry("AC", entry)){
			std::cerr << "AC@" << entry->ID << '\t' << entry->IDX << '\t' << (int)entry->TYPE << std::endl;
		}
		*/


		Iterator::GenotypeIterator it_gt(this->block);
		Iterator::ContainerIteratorDataInterface& temp = *it_gt.iterator_gt_meta.getDataIterator();
		//std::cerr << "Type:" << this->block.gt_support_data_container.header.controller.type << '\t' << this->block.gt_support_data_container.header_stride.controller.type << std::endl;
		//std::cerr << "Size: " << temp.size() << std::endl;
		//std::cerr << this->block.gt_support_data_container.header.stride << '\t' << this->block.gt_support_data_container.header.controller.mixedStride << std::endl;

		// Todo: iterate over GT data
		U32 cost[4];
		cost[0] = 1; cost[1] = 2; cost[2] = 4; cost[3] = 8;
		U64 total_cost = 0;
		for(U32 i = 0; i < temp.size(); ++i){
			const Core::MetaEntry& m = (*it_gt.iterator_meta)[i];
			//std::cerr << it_gt.getCurrentObjectLength() << ',' << it_gt.getCurrentTargetStream() << ' ';
			total_cost += cost[m.hot.controller.gt_primtive_type] * it_gt.getCurrentObjectLength();
			++it_gt;
		}
		//std::cerr << std::endl;
		std::cerr << "Total bytes: " << total_cost << "/" << it_gt.container_rle->getSizeUncompressed() + it_gt.container_simple->getSizeUncompressed() << std::endl;
		assert(total_cost == it_gt.container_rle->getSizeUncompressed() + it_gt.container_simple->getSizeUncompressed());
		it_gt.reset();


		Iterator::MetaIterator* it = this->block.getMetaIterator(); // factory

		// Setup containers
		// INFO
		Iterator::ContainerIterator* info_iterators = nullptr;
		if(this->settings.loadInfoAll){
			info_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_info_streams];

			for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
				info_iterators[i].setup(this->block.info_containers[i]);
			}
		}

		// FORMAT
		Iterator::ContainerIterator* format_iterators = nullptr;
		if(this->settings.loadFormatAll){
			format_iterators = new Iterator::ContainerIterator[this->block.index_entry.n_format_streams];
			for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
				format_iterators[i].setup(this->block.format_containers[i]);
			}
		}

		for(U32 i = 0; i < this->block.size(); ++i){
			const Core::MetaEntry& m = (*it)[i];

			//stream << m.hot.controller.mixed_phasing << ":" << m.hot.controller.phase << ';' << m.hot.controller.anyMissing << ',' << m.hot.controller.anyNA << ':' <<
			//	   m.hot.controller.rle << "->" << m.hot.controller.rle_type << '\t' << (int)m.hot.getGenotypeType() << '\t';

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
					if(this->block.index_entry.filter_bit_vectors[m.getFilterPatternID()][k]){
					// Lookup what that field is
						stream.write(&this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID[0], this->header.getEntry(this->block.index_entry.filter_offsets[k].key).ID.size()) << '\t';
					}
				}
			}

			U32 set = 0;
			if(this->settings.loadInfoAll){
				// Cycle over streams that are set in the given bit-vector

				const Index::BlockIndexBitvector& target_info_vector = this->block.index_entry.info_bit_vectors[m.getInfoPatternID()];
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
				const Index::BlockIndexBitvector& target_format_vector = this->block.index_entry.format_bit_vectors[m.getFormatPatternID()];
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

					set = 0;
					stream.put('\t');
				}

				for(U32 k = 0; k < n_keys_format; ++k){
					if(this->header.entries[this->header.mapTable[this->block.index_entry.format_offsets[firstKey_format[k]].key]].ID == "GT"){
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
		delete [] format_iterators;
		return true;
	}

	U64 iterateGT(std::ostream& stream = std::cout){
		Algorithm::Timer timer;
		timer.Start();
		Iterator::GenotypeIterator it_gt(this->block);

		//Iterator::MetaIterator* it1 = this->block.getMetaIterator();
		//Iterator::MetaIterator* it2 = this->block.getMetaIterator(Core::YON_GT_RLE_DIPLOID_BIALLELIC);
		//Iterator::MetaIterator* it3 = this->block.getMetaIterator(Core::YON_GT_RLE_DIPLOID_NALLELIC);


		//std::cerr << it1->size() << '\t' << it2->size() << '\t' << it3->size() << std::endl;
		//delete it1; delete it2; delete it3;
		//delete it2;

		U32 cost[4];
		cost[0] = 1; cost[1] = 2; cost[2] = 4; cost[3] = 8;
		U64 total_cost = 0;
		for(U32 i = 0; i < this->block.size(); ++i){
			//std::cerr << it_gt.getCurrentObjectLength() << ',' << it_gt.getCurrentTargetStream() << ' ';
			U32 n_sum = 0;
			for(U32 j = 0; j < it_gt.getCurrentObjectLength(); ++j){
				const Core::GTDiploidObject& current_gt = it_gt.getCurrentGTObject();
				//std::cerr << std::bitset<32>(current_gt) << ' ';
				//std::cerr << (int)current_gt.alleles[0] << (current_gt.phase ? '|' : '/') << (int)current_gt.alleles[1] << ':' << current_gt.n_objects << '/' << it_gt.getCurrentObjectLength() << ' ';
				it_gt.incrementGT();
				n_sum += current_gt.n_objects;
			}
			//std::cerr << '\t' << n_sum << std::endl;
			assert(n_sum == 2504);
			total_cost += cost[it_gt.getCurrentMeta().hot.controller.gt_primtive_type] * it_gt.getCurrentObjectLength();
			++it_gt;
		}
		//std::cerr << std::endl;
		//std::cerr << "Total bytes: " << total_cost << "/" << it_gt.container_rle->getSizeUncompressed() + it_gt.container_simple->getSizeUncompressed() << std::endl;
		assert(total_cost == it_gt.container_rle->getSizeUncompressed() + it_gt.container_simple->getSizeUncompressed());
		std::cerr << this->block.size() << '\t' << Helpers::ToPrettyString((U64)((double)this->block.size()*this->header.n_samples/timer.Elapsed().count())) << '\t' << timer.ElapsedString() << std::endl;
		//it_gt.reset();
		return this->block.size();
	}

	bool iterateMeta(std::ostream& stream = std::cout){
		Iterator::MetaIterator* it = this->block.getMetaIterator(); // factory

		//U32 biallelic = 0;
		for(U32 i = 0; i < this->block.size(); ++i){
			//const Core::MetaEntry& m = (*it)[i];
			//biallelic += m.isBiallelic();

		}
		//std::cerr << biallelic << '/' << this->block.size() << std::endl;

		delete it;
		return true;
	}

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	U64 l_data;

	settings_type settings;
	header_type   header;
	index_type    index;
	block_entry_type block;
	codec_manager_type codec_manager;

private:
	U32 n_internal_buffers;
	buffer_type* internal_buffers;
};

}

#endif /* CORE_TACHYONREADER_H_ */
