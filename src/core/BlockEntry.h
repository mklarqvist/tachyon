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
		this->index_entry.reset();
		for(U32 i = 0; i < 100; ++i){
			this->info_containers[i].reset();
			this->format_containers[i].reset();
		}
		this->meta_hot_container.reset();
		this->meta_cold_container.reset();
		this->gt_rle_container.reset();
		this->gt_simple_container.reset();
		this->ppa_manager.reset();

		this->gt_rle_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->gt_simple_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
		this->meta_hot_container.header.controller.type = CORE_TYPE::TYPE_STRUCT;
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
	}

	void updateBaseContainers(buffer_type& buffer){
		this->updateContainer(this->gt_rle_container, this->index_entry.offset_gt_rle, buffer, 0);
		this->updateContainer(this->gt_simple_container, this->index_entry.offset_gt_simple, buffer, 0);
		this->updateContainer(this->meta_hot_container, this->index_entry.offset_hot_meta, buffer, 0);
		this->updateContainer(this->meta_cold_container, this->index_entry.offset_cold_meta, buffer, 0);
	}

	void updateOffsets(void){
		this->index_entry.offset_gt_rle.header = this->gt_rle_container.header;
		this->index_entry.offset_gt_rle.header_stride = this->gt_rle_container.header_stride;
		this->index_entry.offset_gt_simple.header = this->gt_simple_container.header;
		this->index_entry.offset_gt_simple.header_stride = this->gt_simple_container.header_stride;
		this->index_entry.offset_hot_meta.header = this->meta_hot_container.header;
		this->index_entry.offset_cold_meta.header = this->meta_cold_container.header;

		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
			if(this->info_containers[i].header.controller.mixedStride == false) this->index_entry.info_offsets[i].update(this->info_containers[i].header);
			else this->index_entry.info_offsets[i].update(this->info_containers[i].header, this->info_containers[i].header_stride);
		}
		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i){
			if(this->format_containers[i].header.controller.mixedStride == false) this->index_entry.format_offsets[i].update(this->format_containers[i].header);
			else this->index_entry.format_offsets[i].update(this->format_containers[i].header, this->format_containers[i].header_stride);
		}

		}

private:
	void updateContainer(hash_container_type& v, stream_container* container, offset_type* offset, const U32& length, buffer_type& buffer){
		for(U32 i = 0; i < length; ++i){
			if(container[i].buffer_data.size() == 0)
				continue;

			this->updateContainer(container[i], offset[i], buffer, v[i]);
		}
	}

	void updateContainer(stream_container& container, offset_type& offset, buffer_type& buffer, const U32& key){
		offset.update(key);

		if(container.buffer_data.size() == 0)
			return;

		// Check if stream is uniform in content
		if(container.header.controller.type != CORE_TYPE::TYPE_STRUCT)
			container.checkUniformity();

		// Reformat stream to use as small word size as possible
		container.reformat(buffer);

		// Add CRC32 checksum for sanity
		container.generateCRC();

		// Set uncompressed length
		container.header.uLength = container.buffer_data.pointer;

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
		stream << entry.ppa_manager;
		stream << entry.meta_hot_container;
		stream << entry.meta_cold_container;
		stream << entry.gt_rle_container;
		stream << entry.gt_simple_container;

		for(U32 i = 0; i < entry.index_entry.n_info_streams; ++i)
			stream << entry.info_containers[i];

		for(U32 i = 0; i < entry.index_entry.n_format_streams; ++i)
			stream << entry.format_containers[i];

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
