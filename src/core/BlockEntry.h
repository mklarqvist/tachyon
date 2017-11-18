#ifndef CORE_BLOCKENTRY_H_
#define CORE_BLOCKENTRY_H_

namespace Tachyon{
namespace Core{

class BlockEntry{
	typedef BlockEntry self_type;
	typedef Core::StreamContainer stream_container;
	typedef Core::PermutationManager permutation_type;
	typedef Core::EntryHotMetaBase meta_type;
	typedef Core::EntryColdMeta meta_cold_type;
	typedef Index::IndexBlockEntry index_entry_type;

public:
	BlockEntry() :
		info_containers(new stream_container[100]),
		format_containers(new stream_container[100]),
		filter_containers(new stream_container[100])
	{
		for(U32 i = 0; i < 100; ++i){
			this->info_containers[i].resize(65536*4);
			this->format_containers[i].resize(65536*4);
			this->filter_containers[i].resize(65536*4);
		}
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

	// Build Occ function
	// Todo: read and write operations
	void clear(){
		for(U32 i = 0; i < 100; ++i){
			this->info_containers[i].reset();
			this->format_containers[i].reset();
			this->filter_containers[i].reset();
		}
		this->meta_hot_container.reset();
		this->meta_cold_container.reset();
		this->gt_rle_container.reset();
		this->gt_simple_container.reset();
		this->index_entry.reset();
		this->ppa_manager.reset();
	}

	inline void allocateOffsets(const U32& info, const U32& format, const U32& filter){
		//std::cerr << "in allocate: " << this->index_entry.contigID << std::endl;
		//std::cerr << info << '\t' << format << '\t' << filter << std::endl;
		this->index_entry.allocateOffsets(info, format, filter);
	}

	// Sketch:
	// 1) Load index header
	// 2) Load bit vectors
	// 3) Construct local key-value hash table

public:
	index_entry_type  index_entry;
	permutation_type  ppa_manager;
	stream_container  meta_hot_container;
	stream_container  meta_cold_container;
	stream_container  gt_rle_container;
	stream_container  gt_simple_container;
	stream_container* info_containers;
	stream_container* format_containers;
	stream_container* filter_containers;
};

}
}



#endif /* CORE_BLOCKENTRY_H_ */
