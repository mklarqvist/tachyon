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

bool BlockEntry::read(std::ifstream& stream, settings_type& settings){
	settings.load_info_ID_loaded.clear();

	const U64 start_offset = (U64)stream.tellg();
	stream >> this->index_entry;
	const U64 end_of_block = start_offset + this->index_entry.offset_end_of_block;

	if(settings.loadPPA){
		if(this->index_entry.controller.hasGTPermuted){
			stream.seekg(start_offset + this->index_entry.offset_ppa.offset);
			stream >> this->ppa_manager;
		}
	}

	if(settings.loadMetaHot){
		stream.seekg(start_offset + this->index_entry.offset_hot_meta.offset);
		stream >> this->meta_hot_container;
	}

	if(settings.loadMetaCold){
		stream.seekg(start_offset + this->index_entry.offset_cold_meta.offset);
		stream >> this->meta_cold_container;
	}

	if(settings.loadGT){
		stream.seekg(start_offset + this->index_entry.offset_gt_rle.offset);
		stream >> this->gt_rle_container;
	}

	if(settings.loadGTSimple){
		stream.seekg(start_offset + this->index_entry.offset_gt_simple.offset);
		stream >> this->gt_simple_container;
	}

	if(settings.loadInfoAll){
		stream.seekg(start_offset + this->index_entry.info_offsets[0].offset);
		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i)
			stream >> this->info_containers[i];

	} else if(settings.load_info_ID.size() > 0) {
		BlockEntrySettingsMap map_entry;
		// Cycle over all the keys we are interested in
		U32 target_stream_id = 0;
		for(U32 i = 0; i < settings.load_info_ID.size(); ++i){
			// Cycle over all available streams in this block
			for(U32 j = 0; j < this->index_entry.n_info_streams; ++j){
				//std::cerr << "checking: " << this->index_entry.info_offsets[j].key << std::endl;
				if(this->index_entry.info_offsets[j].key == settings.load_info_ID[i]){
					//std::cerr << "Matched: " << i << ',' << j << '\t' << this->index_entry.info_offsets[j].key << "==" << settings.load_info_ID[i] << std::endl;
					settings.load_info_ID_loaded.push_back( BlockEntrySettingsMap(settings.load_info_ID[i], target_stream_id++, j, this->index_entry.info_offsets[j].offset));
					break;
				}
			}
		}

		// Ascertain that random access is linearly forward
		std::sort(settings.load_info_ID_loaded.begin(), settings.load_info_ID_loaded.end());

		// Todo: have to jump to next info block we know exists
		for(U32 i = 0; i < settings.load_info_ID_loaded.size(); ++i){
			//std::cerr << "key: " << settings.load_info_ID_loaded[i].target_stream_local << std::endl;
			//std::cerr << this->index_entry.info_offsets[settings.load_info_ID_loaded[i].target_stream_local].key << '\t' << start_offset + this->index_entry.info_offsets[settings.load_info_ID_loaded[i].target_stream_local].offset << '=' << start_offset + settings.load_info_ID_loaded[i].offset << std::endl;
			stream.seekg(start_offset + settings.load_info_ID_loaded[i].offset);
			if(!stream.good()){
				std::cerr << "failed seek!" << std::endl;
				return false;
			}
			stream >> this->info_containers[i];
			// Store which keys existed
			//std::cerr << "pushing back: " << settings.load_info_ID_loaded[i].target_stream << std::endl;
		}
	} // end case load_info_ID

	if(settings.loadFormatAll){
		stream.seekg(start_offset + this->index_entry.format_offsets[0].offset);
		for(U32 i = 0; i < this->index_entry.n_format_streams; ++i)
			stream >> this->format_containers[i];
	}

	stream.seekg(end_of_block - sizeof(U64));
	U64 eof_marker;
	stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
	assert(eof_marker == Constants::TACHYON_BLOCK_EOF);

	return(true);
}


}
}
