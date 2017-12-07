#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"
#include "BlockEntry.h"

namespace Tachyon{
namespace Core{

BlockEntry::BlockEntry() :
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

BlockEntry::~BlockEntry(){
	delete [] this->info_containers;
}

void BlockEntry::clear(){
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

void BlockEntry::resize(const U32 s){
	if(s == 0) return;
	this->meta_hot_container.resize(s);
	this->meta_cold_container.resize(s);
	this->gt_rle_container.resize(s);
	this->gt_simple_container.resize(s);
}

void BlockEntry::updateBaseContainers(buffer_type& buffer){
	this->updateContainer(this->gt_rle_container,    this->index_entry.offset_gt_rle,    buffer, 0);
	this->updateContainer(this->gt_simple_container, this->index_entry.offset_gt_simple, buffer, 0);
	this->updateContainer(this->meta_hot_container,  this->index_entry.offset_hot_meta,  buffer, 0);
	this->updateContainer(this->meta_cold_container, this->index_entry.offset_cold_meta, buffer, 0);
}

void BlockEntry::updateOffsets(void){
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

void BlockEntry::BlockEntry::updateContainer(hash_container_type& v, stream_container* container, offset_minimal_type* offset, const U32& length, buffer_type& buffer){
	for(U32 i = 0; i < length; ++i){
		offset[i].key = v[i];

		if(container[i].n_entries > 0 && container[i].buffer_data.pointer == 0){
			container[i].header.controller.type = CORE_TYPE::TYPE_BOOLEAN;
			container[i].header.controller.uniform = true;
			container[i].header.stride = 0;
			container[i].header.uLength = 0;
			container[i].header.cLength = 0;
			container[i].header.controller.mixedStride = false;
			container[i].header.controller.encoder = Core::ENCODE_NONE;
			container[i].n_entries = 0;
			container[i].n_additions = 0;
			container[i].header.controller.signedness = 0;
			continue;
		}

		this->updateContainer(container[i], offset[i], buffer, v[i]);
		assert(container[i].header.stride != 0);
	}
}

void BlockEntry::updateContainer(stream_container& container, offset_minimal_type& offset, buffer_type& buffer, const U32& key){
	if(container.buffer_data.size() == 0 && container.header.controller.type != Core::TYPE_BOOLEAN)
		return;

	// Check if stream is uniform in content
	if(container.header.controller.type != CORE_TYPE::TYPE_STRUCT)
		container.checkUniformity();

	// Reformat stream to use as small word size as possible
	container.reformat(buffer);

	// Set uncompressed length
	container.header.uLength = container.buffer_data.pointer;

	// Add CRC32 checksum for sanity
	if(!container.header.controller.mixedStride)
		container.generateCRC(false);

	// If we have mixed striding
	if(container.header.controller.mixedStride){
		// Reformat stream to use as small word size as possible
		container.reformatStride(buffer);
		container.header_stride.uLength = container.buffer_strides.pointer;
		container.generateCRC(true);
	}
}

bool BlockEntry::read(std::ifstream& stream, const settings_type& settings){
	stream >> this->index_entry;
	//std::cerr << this->index_this->contigID << ":" << this->index_this->minPosition << "-" << this->index_this->maxPosition << std::endl;
	//std::cerr << this->index_this->offset_end_of_block - this->index_this->offset_ppa.offset << std::endl;
	//const U64 end_of_block = (U64)stream.tellg() + (this->index_this->offset_end_of_block - this->index_this->offset_ppa.offset) - sizeof(U64);
	//stream.seekg(curpos + (this->index_this->offset_end_of_block - this->index_this->offset_ppa.offset) - sizeof(U64));

	// i.e. search to meta hot
	//const U32 hot_offset = this->index_this->offset_hot_meta.offset - this->index_this->offset_ppa.offset;
	//std::cerr << "seeking to: " << ((U64)stream.tellg() + hot_offset) << std::endl;
	//stream.seekg((U64)stream.tellg() + hot_offset);
	//stream >> this->meta_hot_container;
	//std::cerr << (int)this->meta_hot_container.header.extra[0] << std::endl;
	//std::cerr << "meta data: " << this->meta_hot_container.buffer_data.pointer << std::endl;
	//stream >> this->meta_cold_container;
	//stream >> this->gt_rle_container;
	//std::cerr << this->gt_rle_container.buffer_data.pointer << '\t' << this->index_this->n_variants << std::endl;
	//stream >> this->gt_simple_container;

	//stream.seekg(end_of_block);


	if(settings.loadPPA){
		if(this->index_entry.controller.hasGTPermuted){
			stream >> this->ppa_manager;
		}
	}
	if(settings.loadMetaHot)  stream >> this->meta_hot_container;
	if(settings.loadMetaCold) stream >> this->meta_cold_container;
	if(settings.loadGT)       stream >> this->gt_rle_container;
	if(settings.loadGTSimple) stream >> this->gt_simple_container;
	//stream.seekg(end_of_block);


	//std::cerr << "info_streams: " << this->index_this->n_info_streams << std::endl;
	for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
		stream >> this->info_containers[i];
		//if(this->info_containers[i].header.controller.encoder == ENCODE_NONE){
		//	std::cerr << "ENCODE_NONE | CRC check: " << (this->info_containers[i].checkCRC() ? "PASS" : "FAIL") << std::endl;
		//}
	}

	//std::cerr << "info_streams: " << this->index_this->n_format_streams << std::endl;
	for(U32 i = 0; i < this->index_entry.n_format_streams; ++i)
		stream >> this->format_containers[i];


	U64 eof_marker;
	stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
	assert(eof_marker == Constants::TACHYON_BLOCK_EOF);

	return(true);
}


}
}
