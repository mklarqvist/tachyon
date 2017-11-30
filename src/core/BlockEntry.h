#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

#include "PermutationManager.h"
#include "../index/IndexBlockEntry.h"
#include "StreamContainer.h"
#include "../io/vcf/VCFHeader.h"

namespace Tachyon{
namespace Core{

class BlockEntry{
	typedef BlockEntry self_type;
	typedef Core::StreamContainer stream_container;
	typedef Core::PermutationManager permutation_type;
	typedef Index::IndexBlockEntry index_entry_type;
	typedef Core::Support::HashContainer hash_container_type;
	typedef Core::Support::HashVectorContainer hash_vector_container_type;
	typedef Index::IndexBlockEntryOffsets offset_type;
	typedef Index::IndexBlockEntryHeaderOffsets offset_minimal_type;
	typedef IO::BasicBuffer buffer_type;

public:
	BlockEntry() :
		info_containers(new stream_container[100]),
		format_containers(new stream_container[100])
	{
		for(U32 i = 0; i < 100; ++i){
			this->info_containers[i].resize(65536*4);
			this->format_containers[i].resize(65536*4);
		}

		// Always of type struct
		this->gt_rle_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->gt_simple_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->meta_hot_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->meta_cold_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
	}

	~BlockEntry(){
		delete [] this->info_containers;
	}

	void resize(const U32 s){
		if(s == 0) return;
		this->meta_hot_container.resize(s);
		this->meta_cold_container.resize(s);
		this->gt_rle_container.resize(s);
		this->gt_simple_container.resize(s);
	}

	void clear(){
		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i)
			this->info_containers[i].reset();

		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i)
			this->format_containers[i].reset();

		this->index_entry.reset();
		this->meta_hot_container.reset();
		this->meta_cold_container.reset();
		this->gt_rle_container.reset();
		this->gt_simple_container.reset();
		this->ppa_manager.reset();

		this->gt_rle_container.header.controller.type    = CORE_TYPE::TYPE_STRUCT;
		this->gt_simple_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->meta_hot_container.header.controller.type  = CORE_TYPE::TYPE_STRUCT;
		this->meta_cold_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
	}

	inline void allocateOffsets(const U32& info, const U32& format, const U32& filter){
		this->index_entry.allocateOffsets(info, format, filter);
	}


	inline void updateContainerSet(Index::IndexBlockEntry::INDEX_BLOCK_TARGET target, hash_container_type& v, buffer_type& buffer){
		// Determine target
		switch(target){
		case(Index::IndexBlockEntry::INDEX_BLOCK_TARGET::INDEX_INFO)   :
			return(this->updateContainer(v, this->info_containers, this->index_entry.info_offsets, this->index_entry.n_info_streams, buffer));
			break;
		case(Index::IndexBlockEntry::INDEX_BLOCK_TARGET::INDEX_FORMAT) :
			return(this->updateContainer(v, this->format_containers, this->index_entry.format_offsets, this->index_entry.n_format_streams, buffer));
			break;
		default: std::cerr << "unknown target type" << std::endl; exit(1);
		}
	}

	// For debugging
	void write(std::ofstream& stream, const VCF::VCFHeader& header){
		U64 startPos = stream.tellp();
		stream << this->index_entry;
		std::cout << "Index\t" << (U64)stream.tellp() - startPos << '\n';
		startPos = stream.tellp();
		stream << this->ppa_manager;
		std::cerr << "PPA-v: " << this->ppa_manager.c_length << '\t' << this->ppa_manager.crc << std::endl;
		std::cout << "PPA\t" << (U64)stream.tellp() - startPos << '\n';
		startPos = stream.tellp();
		stream << this->meta_hot_container;
		std::cout << "META-HOT\t" << (U64)stream.tellp() - startPos << '\n';
		startPos = stream.tellp();
		stream << this->meta_cold_container;
		std::cout << "META-COLD\t" << (U64)stream.tellp() - startPos << '\n';
		startPos = stream.tellp();
		stream << this->gt_rle_container;
		std::cout << "GT-RLE\t" << (U64)stream.tellp() - startPos << '\n';
		startPos = stream.tellp();
		stream << this->gt_simple_container;
		std::cout << "GT-SIMPLE\t" << (U64)stream.tellp() - startPos << '\n';

		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
			startPos = stream.tellp();
			stream << this->info_containers[i];
			std::cout << "INFO-" << header[this->index_entry.info_offsets[i].key].ID << '\t' << (U64)stream.tellp() - startPos << '\n';
		}

		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i){
			startPos = stream.tellp();
			stream << this->format_containers[i];
			std::cout << "FORMAT-" << header[this->index_entry.format_offsets[i].key].ID << '\t' << (U64)stream.tellp() - startPos << '\n';
		}

		std::cerr << Helpers::timestamp("DEBUG") << std::endl;
		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i)
			std::cerr << this->index_entry.info_offsets[i].key << '\t';
		std::cerr << std::endl;
		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i)
				std::cerr << this->index_entry.format_offsets[i].key << '\t';
			std::cerr << std::endl;
	}

	void updateBaseContainers(buffer_type& buffer){
		this->updateContainer(this->gt_rle_container, this->index_entry.offset_gt_rle, buffer, 0);
		this->updateContainer(this->gt_simple_container, this->index_entry.offset_gt_simple, buffer, 0);
		this->updateContainer(this->meta_hot_container, this->index_entry.offset_hot_meta, buffer, 0);
		this->updateContainer(this->meta_cold_container, this->index_entry.offset_cold_meta, buffer, 0);
	}

	void updateOffsets(void){
		U32 cum_size = this->index_entry.getDiskSize();
		this->index_entry.offset_ppa.offset = cum_size;
		cum_size += this->ppa_manager.getDiskSize();
		this->index_entry.offset_hot_meta.offset = cum_size;
		cum_size += this->meta_hot_container.getDiskSize();
		this->index_entry.offset_cold_meta.offset = cum_size;
		cum_size += this->meta_cold_container.getDiskSize();
		this->index_entry.offset_gt_rle.offset = cum_size;
		cum_size += this->gt_rle_container.getDiskSize();
		this->index_entry.offset_gt_simple.offset = cum_size;
		cum_size += this->gt_simple_container.getDiskSize();

		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
			this->index_entry.info_offsets[i].offset = cum_size;
			cum_size += this->info_containers[i].getDiskSize();
		}

		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i){
			this->index_entry.format_offsets[i].offset = cum_size;
			cum_size += this->format_containers[i].getDiskSize();
		}

		cum_size += sizeof(U64);
		this->index_entry.offset_end_of_block = cum_size;
		std::cerr << cum_size << std::endl;
	}

private:
	void updateContainer(hash_container_type& v, stream_container* container, offset_minimal_type* offset, const U32& length, buffer_type& buffer){
		for(U32 i = 0; i < length; ++i){
			if(container[i].buffer_data.size() == 0)
				continue;

			this->updateContainer(container[i], offset[i], buffer, v[i]);
		}
	}

	void updateContainer(stream_container& container, offset_minimal_type& offset, buffer_type& buffer, const U32& key){
		offset.key = key;

		if(container.buffer_data.size() == 0)
			return;

		// Check if stream is uniform in content
		if(container.header.controller.type != CORE_TYPE::TYPE_STRUCT)
			container.checkUniformity();

		// Reformat stream to use as small word size as possible
		container.reformat(buffer);

		// Set uncompressed length
		container.header.uLength = container.buffer_data.pointer;

		// Add CRC32 checksum for sanity
		container.generateCRC(true);

		//std::cerr << "update: " << key << '\t' << container.header.uLength << '/' << container.buffer_data.pointer << '\t' << container.header.controller.uniform << std::endl;
		//std::cerr << key << '\t' << container.buffer_data.size() << '\t' << container.header.stride << '\t' << container.header.controller.mixedStride << std::endl;

		// If we have mixed striding
		if(container.header.controller.mixedStride){
			// Reformat stream to use as small word size as possible
			container.reformatStride(buffer);

			container.header_stride.uLength = container.buffer_strides.pointer;
			//std::cerr << key << "-ADD\t" << container.buffer_strides.size() << '\t' << 0 << std::endl;
		}
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.index_entry;
		if(entry.index_entry.controller.hasGTPermuted)
			stream << entry.ppa_manager;
		stream << entry.meta_hot_container;
		stream << entry.meta_cold_container;
		stream << entry.gt_rle_container;
		stream << entry.gt_simple_container;

		for(U32 i = 0; i < entry.index_entry.n_info_streams; ++i)
			stream << entry.info_containers[i];

		for(U32 i = 0; i < entry.index_entry.n_format_streams; ++i)
			stream << entry.format_containers[i];

		stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_BLOCK_EOF), sizeof(U64));

		return(stream);
	}

	// Todo: Temporary!
	// We should only read the entries we're actually interested in!!!!!
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.index_entry;
		//std::cerr << entry.index_entry.contigID << ":" << entry.index_entry.minPosition << "-" << entry.index_entry.maxPosition << std::endl;
		//std::cerr << entry.index_entry.offset_end_of_block - entry.index_entry.offset_ppa.offset << std::endl;
		//const U64 end_of_block = (U64)stream.tellg() + (entry.index_entry.offset_end_of_block - entry.index_entry.offset_ppa.offset) - sizeof(U64);
		//stream.seekg(curpos + (entry.index_entry.offset_end_of_block - entry.index_entry.offset_ppa.offset) - sizeof(U64));

		// i.e. search to meta hot
		//const U32 hot_offset = entry.index_entry.offset_hot_meta.offset - entry.index_entry.offset_ppa.offset;
		//std::cerr << "seeking to: " << ((U64)stream.tellg() + hot_offset) << std::endl;
		//stream.seekg((U64)stream.tellg() + hot_offset);
		//stream >> entry.meta_hot_container;
		//std::cerr << (int)entry.meta_hot_container.header.extra[0] << std::endl;
		//std::cerr << "meta data: " << entry.meta_hot_container.buffer_data.pointer << std::endl;
		//stream >> entry.meta_cold_container;
		//stream >> entry.gt_rle_container;
		//std::cerr << entry.gt_rle_container.buffer_data.pointer << '\t' << entry.index_entry.n_variants << std::endl;
		//stream >> entry.gt_simple_container;

		//stream.seekg(end_of_block);


		if(entry.index_entry.controller.hasGTPermuted) stream >> entry.ppa_manager;
		stream >> entry.meta_hot_container;
		stream >> entry.meta_cold_container;
		stream >> entry.gt_rle_container;
		stream >> entry.gt_simple_container;
		//stream.seekg(end_of_block);


		//std::cerr << "info_streams: " << entry.index_entry.n_info_streams << std::endl;
		for(U32 i = 0; i < entry.index_entry.n_info_streams; ++i){
			stream >> entry.info_containers[i];
			//if(entry.info_containers[i].header.controller.encoder == ENCODE_NONE){
			//	std::cerr << "ENCODE_NONE | CRC check: " << (entry.info_containers[i].checkCRC() ? "PASS" : "FAIL") << std::endl;
			//}
		}

		//std::cerr << "info_streams: " << entry.index_entry.n_format_streams << std::endl;
		for(U32 i = 0; i < entry.index_entry.n_format_streams; ++i)
			stream >> entry.format_containers[i];


		U64 eof_marker;
		stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
		assert(eof_marker == Constants::TACHYON_BLOCK_EOF);


		return(stream);
	}

public:
	index_entry_type  index_entry;
	permutation_type  ppa_manager;
	stream_container  meta_hot_container;
	stream_container  meta_cold_container;
	stream_container  gt_rle_container;
	stream_container  gt_simple_container;
	stream_container* info_containers;
	stream_container* format_containers;
};

}
}

#endif /* CORE_BLOCKENTRY_H_ */
