#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"
#include "Block.h"

namespace Tachyon{
namespace Core{

Block::Block() :
	info_containers(new container_type[200]),
	format_containers(new container_type[200])
{
	// Base container streams are always of type TYPE_STRUCT
	this->meta_format_map_ids.setType(YON_TYPE_32B);
	this->meta_filter_map_ids.setType(YON_TYPE_32B);
	this->meta_info_map_ids.setType(YON_TYPE_32B);
	this->gt_rle_container.setType(YON_TYPE_STRUCT);
	this->gt_simple_container.setType(YON_TYPE_STRUCT);
	this->gt_support_data_container.setType(YON_TYPE_32B);
	this->gt_support_data_container.header.controller.signedness = true;
	this->meta_hot_container.setType(YON_TYPE_STRUCT);
	this->meta_cold_container.setType(YON_TYPE_STRUCT);
	this->meta_info_map_ids.setStrideSize(1);
	this->meta_filter_map_ids.setStrideSize(1);
	this->meta_format_map_ids.setStrideSize(1);
	this->gt_support_data_container.setStrideSize(1);
}

Block::~Block(){
	delete [] this->info_containers;
}

void Block::clear(void){
	for(U32 i = 0; i < this->index_entry.n_info_streams; ++i)
		this->info_containers[i].reset();

	for(U32 i = 0; i < this->index_entry.n_format_streams; ++i)
		this->format_containers[i].reset();

	this->index_entry.reset();
	this->meta_info_map_ids.reset();
	this->meta_filter_map_ids.reset();
	this->meta_format_map_ids.reset();
	this->meta_hot_container.reset();
	this->meta_cold_container.reset();
	this->gt_support_data_container.reset();
	this->gt_rle_container.reset();
	this->gt_simple_container.reset();
	this->ppa_manager.reset();

	// Base container data types are always TYPE_STRUCT
	// Map ID fields are always U32 fields
	this->meta_format_map_ids.setType(YON_TYPE_32B);
	this->meta_filter_map_ids.setType(YON_TYPE_32B);
	this->meta_info_map_ids.setType(YON_TYPE_32B);
	this->gt_rle_container.setType(YON_TYPE_STRUCT);
	this->gt_simple_container.setType(YON_TYPE_STRUCT);
	this->gt_support_data_container.setType(YON_TYPE_32B);
	this->meta_hot_container.setType(YON_TYPE_STRUCT);
	this->meta_cold_container.setType(YON_TYPE_STRUCT);
	this->meta_info_map_ids.setStrideSize(1);
	this->meta_filter_map_ids.setStrideSize(1);
	this->meta_format_map_ids.setStrideSize(1);
	this->gt_support_data_container.setStrideSize(1);
}

void Block::resize(const U32 s){
	if(s == 0) return;
	this->meta_hot_container.resize(s);
	this->meta_cold_container.resize(s);
	this->gt_rle_container.resize(s);
	this->gt_simple_container.resize(s);
	this->meta_info_map_ids.resize(s);
	this->meta_filter_map_ids.resize(s);
	this->meta_format_map_ids.resize(s);
	this->gt_support_data_container.resize(s);
	for(U32 i = 0; i < 200; ++i){
		this->info_containers[i].resize(s);
		this->format_containers[i].resize(s);
	}
}

void Block::updateBaseContainers(void){
	this->updateContainer(this->gt_rle_container);
	this->updateContainer(this->gt_simple_container);
	this->updateContainer(this->gt_support_data_container);
	this->updateContainer(this->meta_hot_container);
	this->updateContainer(this->meta_cold_container);
	this->updateContainer(this->meta_filter_map_ids, false);
	this->updateContainer(this->meta_format_map_ids, false);
	this->updateContainer(this->meta_info_map_ids  , false);
}

void Block::updateOffsets(void){
	U32 cum_size = this->index_entry.getObjectSize();
	this->index_entry.offset_ppa.offset = cum_size;
	cum_size += this->ppa_manager.getObjectSize();

	this->index_entry.offset_hot_meta.offset = cum_size;
	cum_size += this->meta_hot_container.getObjectSize();

	this->index_entry.offset_cold_meta.offset = cum_size;
	cum_size += this->meta_cold_container.getObjectSize();

	this->index_entry.offset_gt_rle.offset = cum_size;
	cum_size += this->gt_rle_container.getObjectSize();

	this->index_entry.offset_gt_simple.offset = cum_size;
	cum_size += this->gt_simple_container.getObjectSize();


	this->index_entry.offset_gt_helper.offset = cum_size;
	cum_size += this->gt_support_data_container.getObjectSize();

	this->index_entry.offset_meta_info_id.offset = cum_size;
	cum_size += this->meta_info_map_ids.getObjectSize();

	this->index_entry.offset_meta_filter_id.offset = cum_size;
	cum_size += this->meta_filter_map_ids.getObjectSize();

	this->index_entry.offset_meta_format_id.offset = cum_size;
	cum_size += this->meta_format_map_ids.getObjectSize();


	for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
		this->index_entry.info_offsets[i].offset = cum_size;
		cum_size += this->info_containers[i].getObjectSize();
	}

	for(U32 i = 0; i < this->index_entry.n_format_streams; ++i){
		this->index_entry.format_offsets[i].offset = cum_size;
		cum_size += this->format_containers[i].getObjectSize();
	}

	// Size of EOF marker
	cum_size += sizeof(U64);

	// Update final size
	this->index_entry.offset_end_of_block = cum_size;
}

void Block::Block::updateContainer(container_type* container, const U32& length){
	for(U32 i = 0; i < length; ++i){
		// If the data container has entries in it but has
		// no actual data then it is a BOOLEAN
		if(container[i].n_entries > 0 && container[i].buffer_data_uncompressed.pointer == 0){
			container[i].header.controller.type = YON_TYPE_BOOLEAN;
			container[i].header.controller.uniform = true;
			container[i].header.stride  = 0;
			container[i].header.uLength = 0;
			container[i].header.cLength = 0;
			container[i].header.controller.mixedStride = false;
			container[i].header.controller.encoder = Core::YON_ENCODE_NONE;
			container[i].n_entries      = 0;
			container[i].n_additions    = 0;
			container[i].header.controller.signedness = 0;
			continue;
		}

		// If the data type is not a BOOLEAN
		this->updateContainer(container[i]);
		assert(container[i].header.stride != 0);
	}
}

void Block::updateContainer(container_type& container, bool reformat){
	if(container.buffer_data_uncompressed.size() == 0 &&
	   container.header.controller.type != Core::YON_TYPE_BOOLEAN)
		return;

	// Check if stream is uniform in content
	if(container.header.controller.type != YON_TYPE_STRUCT){
		container.checkUniformity();
		// Reformat stream to use as small word size as possible
		if(reformat) container.reformat();
	}

	// Set uncompressed length
	container.header.uLength = container.buffer_data_uncompressed.pointer;

	// If we have mixed striding
	if(container.header.controller.mixedStride == true){
		// Reformat stream to use as small word size as possible
		if(reformat) container.reformatStride();
		container.header_stride.uLength = container.buffer_strides_uncompressed.pointer;
	}
}

bool Block::read(std::ifstream& stream, settings_type& settings){
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

	if(settings.loadGTSimple || settings.loadGT){
		stream.seekg(start_offset + this->index_entry.offset_gt_helper.offset);
		stream >> this->gt_support_data_container;
	}

	if(settings.loadInfoAll || settings.load_info_ID_loaded.size()){
		stream.seekg(start_offset + this->index_entry.offset_meta_info_id.offset);
		stream >> this->meta_info_map_ids;
	}

	// Todo: always load?
	stream.seekg(start_offset + this->index_entry.offset_meta_filter_id.offset);
	stream >> this->meta_filter_map_ids;

	if(settings.loadInfoAll){
		stream.seekg(start_offset + this->index_entry.offset_meta_format_id.offset);
		stream >> this->meta_format_map_ids;
	}

	// Load all info
	if(settings.loadInfoAll){
		stream.seekg(start_offset + this->index_entry.info_offsets[0].offset);
		for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
			stream >> this->info_containers[i];
			//std::cerr << "loaded: " << this->index_entry.info_offsets[i].key << '\t' << this->info_containers[i].header.cLength << std::endl;
		}

	}
	// If we have supplied a list of identifiers
	else if(settings.load_info_ID.size() > 0) {
		BlockEntrySettingsMap map_entry;
		// Cycle over all the keys we are interested in
		U32 iterator_index = 0;
		for(U32 i = 0; i < settings.load_info_ID.size(); ++i){
			// Cycle over all available streams in this block
			for(U32 j = 0; j < this->index_entry.n_info_streams; ++j){
				// If there is a match
				// Push back field into map
				if(this->index_entry.info_offsets[j].key == settings.load_info_ID[i]){
					settings.load_info_ID_loaded.push_back(
							BlockEntrySettingsMap(
									iterator_index++,                   // iterator value
									j,                                  // local index id
									&this->index_entry.info_offsets[j]) // offset
									);
					break;
				}
			}
		}

		// Ascertain that random access is linearly forward
		std::sort(settings.load_info_ID_loaded.begin(), settings.load_info_ID_loaded.end());

		// Todo: have to jump to next info block we know exists
		for(U32 i = 0; i < settings.load_info_ID_loaded.size(); ++i){
			stream.seekg(start_offset + settings.load_info_ID_loaded[i].offset->offset);
			if(!stream.good()){
				std::cerr << Helpers::timestamp("ERROR","IO") << "Failed seek!" << std::endl;
				return false;
			}

			// Read data
			stream >> this->info_containers[settings.load_info_ID_loaded[i].iterator_index];
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

bool Block::write(std::ofstream& stream,
                   import_stats_type& stats,
                   import_stats_type& stats_uncompressed){
	U64 last_pos = stream.tellp();
	stream << this->index_entry;
	stats.total_header_cost += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	if(this->index_entry.controller.hasGTPermuted){
		stream << this->ppa_manager;
		stats.total_ppa_cost += (U64)stream.tellp() - last_pos;
		stats_uncompressed.total_ppa_cost += this->ppa_manager.u_length;
		last_pos = stream.tellp();
	}

	stream << this->meta_hot_container;
	stream << this->meta_cold_container;
	stats.total_meta_cost += (U64)stream.tellp() - last_pos;
	stats_uncompressed.total_meta_cost += this->meta_hot_container.buffer_data_uncompressed.pointer;
	stats_uncompressed.total_meta_cost += this->meta_cold_container.buffer_data_uncompressed.pointer;
	last_pos = stream.tellp();

	stream << this->gt_rle_container;
	stream << this->gt_simple_container;
	stream << this->gt_support_data_container;
	stats.total_gt_cost += (U64)stream.tellp() - last_pos;
	stats_uncompressed.total_gt_cost += this->gt_rle_container.buffer_data_uncompressed.pointer;
	stats_uncompressed.total_gt_cost += this->gt_simple_container.buffer_data_uncompressed.pointer;
	stats_uncompressed.total_gt_cost += this->gt_support_data_container.buffer_data_uncompressed.pointer;
	last_pos = stream.tellp();

	stream << this->meta_info_map_ids;
	stream << this->meta_filter_map_ids;
	stream << this->meta_format_map_ids;
	stats.total_special_cost += (U64)stream.tellp() - last_pos;
	stats_uncompressed.total_special_cost += this->meta_info_map_ids.buffer_data_uncompressed.pointer;
	stats_uncompressed.total_special_cost += this->meta_filter_map_ids.buffer_data_uncompressed.pointer;
	stats_uncompressed.total_special_cost += this->meta_format_map_ids.buffer_data_uncompressed.pointer;
	last_pos = stream.tellp();

	for(U32 i = 0; i < this->index_entry.n_info_streams; ++i){
		stream << this->info_containers[i];
		stats_uncompressed.total_info_cost += this->info_containers[i].buffer_data_uncompressed.pointer;
	}

	stats.total_info_cost += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	for(U32 i = 0; i < this->index_entry.n_format_streams; ++i){
		stream << this->format_containers[i];
		stats_uncompressed.total_format_cost += this->format_containers[i].buffer_data_uncompressed.pointer;
	}

	stats.total_format_cost += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	stream.write(reinterpret_cast<const char*>(&Constants::TACHYON_BLOCK_EOF), sizeof(U64));

	return(true);
}

Iterator::MetaIterator* Block::getMetaIterator(void) const{
	meta_iterator_type* it;
	if(this->meta_cold_container.buffer_data_uncompressed.size())
		it = new meta_iterator_type(this->meta_hot_container, this->meta_cold_container);
	else
		it = new meta_iterator_type(this->meta_hot_container);

	if(this->meta_info_map_ids.buffer_data_uncompressed.size())
		it->setInfoIDContainer(this->meta_info_map_ids);

	if(this->meta_filter_map_ids.buffer_data_uncompressed.size())
		it->setFilterIDContainer(this->meta_filter_map_ids);

	if(this->meta_format_map_ids.buffer_data_uncompressed.size())
		it->setFormatIDContainer(this->meta_format_map_ids);

	return(it);
}

Iterator::MetaIterator* Block::getMetaIterator(const Core::TACHYON_GT_TYPE gt_filter) const{
	meta_iterator_type* it;
	if(this->meta_cold_container.buffer_data_uncompressed.size())
		it = new meta_iterator_type(this->meta_hot_container, this->meta_cold_container, gt_filter);
	else
		it = new meta_iterator_type(this->meta_hot_container, gt_filter);

	if(this->meta_info_map_ids.buffer_data_uncompressed.size())
		it->setInfoIDContainer(this->meta_info_map_ids);

	if(this->meta_filter_map_ids.buffer_data_uncompressed.size())
		it->setFilterIDContainer(this->meta_filter_map_ids);

	if(this->meta_format_map_ids.buffer_data_uncompressed.size())
		it->setFormatIDContainer(this->meta_format_map_ids);

	return(it);
}

}
}
