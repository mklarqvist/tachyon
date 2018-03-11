#include "variantblock.h"

#include "../support/helpers.h"
#include "../support/type_definitions.h"

namespace tachyon{
namespace containers{

VariantBlock::VariantBlock() :
	info_containers(new container_type[200]),
	format_containers(new container_type[200])
{
	// Base container streams are always of type TYPE_STRUCT
	this->meta_format_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_format_map_ids.header.data_header.controller.signedness = 1;
	this->meta_filter_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_filter_map_ids.header.data_header.controller.signedness = 1;
	this->meta_info_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_info_map_ids.header.data_header.controller.signedness = 1;
	this->gt_rle_container.setType(tachyon::YON_TYPE_STRUCT);
	this->gt_simple_container.setType(tachyon::YON_TYPE_STRUCT);
	this->gt_support_data_container.setType(tachyon::YON_TYPE_32B);
	this->gt_support_data_container.header.data_header.controller.signedness = true;
	this->meta_hot_container.setType(tachyon::YON_TYPE_STRUCT);
	this->meta_cold_container.setType(tachyon::YON_TYPE_STRUCT);
	this->meta_info_map_ids.setStrideSize(1);
	this->meta_filter_map_ids.setStrideSize(1);
	this->meta_format_map_ids.setStrideSize(1);
	this->gt_support_data_container.setStrideSize(1);
}

VariantBlock::~VariantBlock(){
	delete [] this->info_containers;
	delete [] this->format_containers;
}

void VariantBlock::clear(void){
	for(U32 i = 0; i < this->header.n_info_streams; ++i)
		this->info_containers[i].reset();

	for(U32 i = 0; i < this->header.n_format_streams; ++i)
		this->format_containers[i].reset();

	this->header.reset();
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
	// Map ID fields are always S32 fields
	this->meta_format_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_format_map_ids.header.data_header.controller.signedness = 1;
	this->meta_filter_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_filter_map_ids.header.data_header.controller.signedness = 1;
	this->meta_info_map_ids.setType(tachyon::YON_TYPE_32B);
	this->meta_info_map_ids.header.data_header.controller.signedness = 1;
	this->gt_rle_container.setType(tachyon::YON_TYPE_STRUCT);
	this->gt_simple_container.setType(tachyon::YON_TYPE_STRUCT);
	this->gt_support_data_container.setType(tachyon::YON_TYPE_32B);
	this->meta_hot_container.setType(tachyon::YON_TYPE_STRUCT);
	this->meta_cold_container.setType(tachyon::YON_TYPE_STRUCT);
	this->meta_info_map_ids.setStrideSize(1);
	this->meta_filter_map_ids.setStrideSize(1);
	this->meta_format_map_ids.setStrideSize(1);
	this->gt_support_data_container.setStrideSize(1);
}

void VariantBlock::resize(const U32 s){
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

void VariantBlock::updateBaseContainers(void){
	this->updateContainer(this->gt_rle_container);
	this->updateContainer(this->gt_simple_container);
	this->updateContainer(this->gt_support_data_container);
	this->updateContainer(this->meta_hot_container);
	this->updateContainer(this->meta_cold_container);
	this->updateContainer(this->meta_filter_map_ids, false);
	this->updateContainer(this->meta_format_map_ids, false);
	this->updateContainer(this->meta_info_map_ids  , false);
}

void VariantBlock::updateOffsets(void){
	U32 cum_size = this->header.getObjectSize();

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
		this->header.offset_ppa.data_header.offset = cum_size;
		cum_size += this->ppa_manager.getObjectSize();
	}

	this->header.offset_hot_meta = this->meta_hot_container.header; // move header over
	this->header.offset_hot_meta.data_header.offset = cum_size;
	cum_size += this->meta_hot_container.getObjectSize();

	this->header.offset_cold_meta = this->meta_cold_container.header; // move header over
	this->header.offset_cold_meta.data_header.offset = cum_size;
	cum_size += this->meta_cold_container.getObjectSize();

	if(this->gt_rle_container.buffer_data_uncompressed.size()){
		this->header.offset_gt_rle = this->gt_rle_container.header; // move header over
		this->header.offset_gt_rle.data_header.offset = cum_size;
		cum_size += this->gt_rle_container.getObjectSize();
	}

	if(this->gt_simple_container.buffer_data_uncompressed.size()){
		this->header.offset_gt_simple = this->gt_simple_container.header; // move header over
		this->header.offset_gt_simple.data_header.offset = cum_size;
		cum_size += this->gt_simple_container.getObjectSize();
	}

	if(this->gt_support_data_container.buffer_data_uncompressed.size()){
		this->header.offset_gt_helper = this->gt_support_data_container.header; // move header over
		this->header.offset_gt_helper.data_header.offset = cum_size;
		cum_size += this->gt_support_data_container.getObjectSize();
	}

	if(this->meta_info_map_ids.buffer_data_uncompressed.size()){
		this->header.offset_meta_info_id = this->meta_info_map_ids.header; // move header over
		this->header.offset_meta_info_id.data_header.offset = cum_size;
		cum_size += this->meta_info_map_ids.getObjectSize();
	}

	if(this->meta_filter_map_ids.buffer_data_uncompressed.size()){
		this->header.offset_meta_filter_id = this->meta_filter_map_ids.header; // move header over
		this->header.offset_meta_filter_id.data_header.offset = cum_size;
		cum_size += this->meta_filter_map_ids.getObjectSize();
	}

	if(this->meta_format_map_ids.buffer_data_uncompressed.size()){
		this->header.offset_meta_format_id = this->meta_format_map_ids.header; // move header over
		this->header.offset_meta_format_id.data_header.offset = cum_size;
		cum_size += this->meta_format_map_ids.getObjectSize();
	}

	for(U32 i = 0; i < this->header.n_info_streams; ++i){
		this->header.info_offsets[i] = this->info_containers[i].header; // move header over
		this->header.info_offsets[i].data_header.offset = cum_size;
		cum_size += this->info_containers[i].getObjectSize();
	}

	for(U32 i = 0; i < this->header.n_format_streams; ++i){
		this->header.format_offsets[i] = this->format_containers[i].header; // move header over
		this->header.format_offsets[i].data_header.offset = cum_size;
		cum_size += this->format_containers[i].getObjectSize();
	}

	// Size of EOF marker
	cum_size += sizeof(U64);

	// Update final size
	this->header.offset_end_of_block = cum_size;
}

void VariantBlock::VariantBlock::updateContainer(container_type* container, const U32& length){
	for(U32 i = 0; i < length; ++i){
		// If the data container has entries in it but has
		// no actual data then it is a BOOLEAN
		if(container[i].n_entries > 0 && container[i].buffer_data_uncompressed.size() == 0){
			container[i].header.data_header.controller.type = tachyon::YON_TYPE_BOOLEAN;
			container[i].header.data_header.controller.uniform = true;
			container[i].header.data_header.stride  = 0;
			container[i].header.data_header.uLength = 0;
			container[i].header.data_header.cLength = 0;
			container[i].header.data_header.controller.mixedStride = false;
			container[i].header.data_header.controller.encoder = tachyon::YON_ENCODE_NONE;
			container[i].n_entries      = 0;
			container[i].n_additions    = 0;
			container[i].header.data_header.controller.signedness = 0;
			continue;
		}

		// If the data type is not a BOOLEAN
		this->updateContainer(container[i]);
		assert(container[i].header.data_header.stride != 0);
	}
}

void VariantBlock::updateContainer(container_type& container, bool reformat){
	if(container.buffer_data_uncompressed.size() == 0 &&
	   container.header.data_header.controller.type != tachyon::YON_TYPE_BOOLEAN)
		return;

	// Check if stream is uniform in content
	if(container.header.data_header.controller.type != tachyon::YON_TYPE_STRUCT){
		container.checkUniformity();
		// Reformat stream to use as small word size as possible
		if(reformat) container.reformat();
	}

	// Set uncompressed length
	container.header.data_header.uLength = container.buffer_data_uncompressed.size();

	// If we have mixed striding
	if(container.header.data_header.hasMixedStride()){
		// Reformat stream to use as small word size as possible
		if(reformat) container.reformatStride();
		container.header.stride_header.uLength = container.buffer_strides_uncompressed.size();
	}
}

bool VariantBlock::read(std::ifstream& stream, settings_type& settings){
	settings.load_info_ID_loaded.clear();

	const U64 start_offset = (U64)stream.tellg();
	stream >> this->header;
	const U64 end_of_block = start_offset + this->header.offset_end_of_block;

	if(settings.importPPA){
		if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
			stream.seekg(start_offset + this->header.offset_ppa.data_header.offset);
			stream >> this->ppa_manager;
		}
	}

	if(settings.importMetaHot){
		stream.seekg(start_offset + this->header.offset_hot_meta.data_header.offset);
		stream >> this->meta_hot_container;
		this->meta_hot_container.header = this->header.offset_hot_meta;
	}

	if(settings.importMetaCold){
		stream.seekg(start_offset + this->header.offset_cold_meta.data_header.offset);
		stream >> this->meta_cold_container;
		this->meta_cold_container.header = this->header.offset_cold_meta;
	}

	if(settings.importGT && this->header.offset_gt_rle.data_header.offset != -1){
		stream.seekg(start_offset + this->header.offset_gt_rle.data_header.offset);
		stream >> this->gt_rle_container;
		this->gt_rle_container.header = this->header.offset_gt_rle;
	}

	if(settings.importGTSimple && this->header.offset_gt_simple.data_header.offset != -1){
		stream.seekg(start_offset + this->header.offset_gt_simple.data_header.offset);
		stream >> this->gt_simple_container;
		this->gt_simple_container.header = this->header.offset_gt_simple;
	}

	if((settings.importGTSimple || settings.importGT) && this->header.offset_gt_helper.data_header.offset != -1){
		stream.seekg(start_offset + this->header.offset_gt_helper.data_header.offset);
		stream >> this->gt_support_data_container;
	}

	if((settings.importMetaHot || settings.importMetaCold) && this->header.offset_meta_info_id.data_header.offset != -1){
		stream.seekg(start_offset + this->header.offset_meta_info_id.data_header.offset);
		stream >> this->meta_info_map_ids;
		this->meta_info_map_ids.header = this->header.offset_meta_info_id;
		stream >> this->meta_filter_map_ids;
		this->meta_filter_map_ids.header = this->header.offset_meta_filter_id;
		stream >> this->meta_format_map_ids;
		this->meta_format_map_ids.header = this->header.offset_meta_format_id;
	}

	// Load all info
	if(settings.importInfoAll && this->header.n_info_streams > 0){
		stream.seekg(start_offset + this->header.info_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->header.n_info_streams; ++i){
			stream >> this->info_containers[i];
			this->info_containers[i].header = this->header.info_offsets[i];
			//std::cerr << "loaded: " << this->index_entry.info_offsets[i].global_key << '\t' << this->info_containers[i].header.cLength << std::endl;
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.info_ID_list.size() > 0) {
		core::SettingsMap map_entry;
		// Cycle over all the keys we are interested in
		U32 iterator_index = 0;
		for(U32 i = 0; i < settings.info_ID_list.size(); ++i){
			// Cycle over all available streams in this block
			for(U32 j = 0; j < this->header.n_info_streams; ++j){
				// If there is a match
				// Push back field into map
				if(this->header.info_offsets[j].data_header.global_key == settings.info_ID_list[i]){
					settings.load_info_ID_loaded.push_back(
							core::SettingsMap(
									iterator_index++,                   // iterator value
									j,                                  // local index id
									&this->header.info_offsets[j]) // offset
									);
					break;
				}
			}
		}

		// Ascertain that random access is linearly forward
		std::sort(settings.load_info_ID_loaded.begin(), settings.load_info_ID_loaded.end());

		// Todo: have to jump to next info block we know exists
		for(U32 i = 0; i < settings.load_info_ID_loaded.size(); ++i){
			stream.seekg(start_offset + settings.load_info_ID_loaded[i].offset->data_header.offset);
			if(!stream.good()){
				std::cerr << utility::timestamp("ERROR","IO") << "Failed seek!" << std::endl;
				return false;
			}

			// Read data
			stream >> this->info_containers[settings.load_info_ID_loaded[i].iterator_index];
			this->info_containers[settings.load_info_ID_loaded[i].iterator_index].header = this->header.info_offsets[settings.load_info_ID_loaded[i].iterator_index];
		}
	} // end case load_info_ID

	if(settings.importFormatAll && this->header.n_format_streams){
		stream.seekg(start_offset + this->header.format_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->header.n_format_streams; ++i){
			stream >> this->format_containers[i];
			this->format_containers[i].header = this->header.format_offsets[i];
			//std::cerr << "loaded: " << this->index_entry.format_offsets[i].global_key << '\t' << this->format_containers[i].header.cLength << std::endl;
		}
	}

	stream.seekg(end_of_block - sizeof(U64));
	U64 eof_marker;
	stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
	assert(eof_marker == constants::TACHYON_BLOCK_EOF);
	return(true);
}

bool VariantBlock::write(std::ofstream& stream,
                         import_stats_type& stats_basic,
                         import_stats_type& stats_info,
                         import_stats_type& stats_format)
{
	U64 last_pos = stream.tellp();
	const U64 start_pos = last_pos;
	stream << this->header;
	stats_basic[0].cost_uncompressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
		stream << this->ppa_manager;
		stats_basic[1].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[1].cost_uncompressed += this->ppa_manager.getObjectSize();
		last_pos = stream.tellp();
	}

	stream << this->meta_hot_container;
	stats_basic[2].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[2].cost_uncompressed += this->meta_hot_container.getObjectSize();
	stats_basic[2].cost_uncompressed += (S32)this->meta_hot_container.header.data_header.uLength - this->meta_hot_container.header.data_header.cLength;
	stats_basic[2].cost_uncompressed += (S32)this->meta_hot_container.header.stride_header.uLength - this->meta_hot_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	stream << this->meta_cold_container;
	stats_basic[3].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[3].cost_uncompressed += this->meta_cold_container.getObjectSize();
	stats_basic[3].cost_uncompressed += (S32)this->meta_cold_container.header.data_header.uLength - this->meta_cold_container.header.data_header.cLength;
	stats_basic[3].cost_uncompressed += (S32)this->meta_cold_container.header.stride_header.uLength - this->meta_cold_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	if(this->gt_rle_container.getSizeUncompressed()){
		stream << this->gt_rle_container;
		stats_basic[4].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[4].cost_uncompressed += this->gt_rle_container.getObjectSize();
		stats_basic[4].cost_uncompressed += (S32)this->gt_rle_container.header.data_header.uLength - this->gt_rle_container.header.data_header.cLength;
		stats_basic[4].cost_uncompressed += (S32)this->gt_rle_container.header.stride_header.uLength - this->gt_rle_container.header.stride_header.cLength;
		last_pos = stream.tellp();
	}

	if(this->gt_simple_container.getSizeUncompressed()){
		stream << this->gt_simple_container;
		stats_basic[5].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[5].cost_uncompressed += this->gt_simple_container.getObjectSize();
		stats_basic[5].cost_uncompressed += (S32)this->gt_simple_container.header.data_header.uLength - this->gt_simple_container.header.data_header.cLength;
		stats_basic[5].cost_uncompressed += (S32)this->gt_simple_container.header.stride_header.uLength - this->gt_simple_container.header.stride_header.cLength;
		last_pos = stream.tellp();
	}

	if(this->gt_support_data_container.getSizeUncompressed()){
		stream << this->gt_support_data_container;
		stats_basic[6].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[6].cost_uncompressed += this->gt_support_data_container.getObjectSize();
		stats_basic[6].cost_uncompressed += (S32)this->gt_support_data_container.header.data_header.uLength - this->gt_support_data_container.header.data_header.cLength;
		stats_basic[6].cost_uncompressed += (S32)this->gt_support_data_container.header.stride_header.uLength - this->gt_support_data_container.header.stride_header.cLength;
		last_pos = stream.tellp();
	}

	if(this->meta_info_map_ids.getSizeUncompressed()){
		stream << this->meta_info_map_ids;
		stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[7].cost_uncompressed += this->meta_info_map_ids.getObjectSize();
		stats_basic[7].cost_uncompressed += (S32)this->meta_info_map_ids.header.data_header.uLength - this->meta_info_map_ids.header.data_header.cLength;
		stats_basic[7].cost_uncompressed += (S32)this->meta_info_map_ids.header.stride_header.uLength - this->meta_info_map_ids.header.stride_header.cLength;

		last_pos = stream.tellp();
	}

	if(this->meta_filter_map_ids.getSizeUncompressed()){
		stream << this->meta_filter_map_ids;
		stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[7].cost_uncompressed += this->meta_filter_map_ids.getObjectSize();
		stats_basic[7].cost_uncompressed += (S32)this->meta_filter_map_ids.header.data_header.uLength - this->meta_filter_map_ids.header.data_header.cLength;
		stats_basic[7].cost_uncompressed += (S32)this->meta_filter_map_ids.header.stride_header.uLength - this->meta_filter_map_ids.header.stride_header.cLength;

		last_pos = stream.tellp();
	}

	if(this->meta_format_map_ids.getSizeUncompressed()){
		stream << this->meta_format_map_ids;
		stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[7].cost_uncompressed += this->meta_format_map_ids.getObjectSize();
		stats_basic[7].cost_uncompressed += (S32)this->meta_format_map_ids.header.data_header.uLength - this->meta_format_map_ids.header.data_header.cLength;
		stats_basic[7].cost_uncompressed += (S32)this->meta_format_map_ids.header.stride_header.uLength - this->meta_format_map_ids.header.stride_header.cLength;

		last_pos = stream.tellp();
	}

	for(U32 i = 0; i < this->header.n_info_streams; ++i){
		stream << this->info_containers[i];
		stats_info[this->header.info_offsets[i].data_header.global_key].cost_uncompressed += this->info_containers[i].header.data_header.uLength;
		stats_info[this->header.info_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->info_containers[i].header.data_header.uLength - this->info_containers[i].header.data_header.cLength;
		stats_info[this->header.info_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->info_containers[i].header.stride_header.uLength - this->info_containers[i].header.stride_header.cLength;
		stats_info[this->header.info_offsets[i].data_header.global_key].cost_compressed   += this->info_containers[i].header.data_header.cLength;
		stats_basic[8].cost_uncompressed += this->info_containers[i].getObjectSize();
		stats_basic[8].cost_uncompressed += (S32)this->info_containers[i].header.data_header.uLength - this->info_containers[i].header.data_header.cLength;
		stats_basic[8].cost_uncompressed += (S32)this->info_containers[i].header.stride_header.uLength - this->info_containers[i].header.stride_header.cLength;
	}

	stats_basic[8].cost_compressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	for(U32 i = 0; i < this->header.n_format_streams; ++i){
		stream << this->format_containers[i];
		stats_format[this->header.format_offsets[i].data_header.global_key].cost_uncompressed += this->format_containers[i].header.data_header.uLength;
		stats_format[this->header.format_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->format_containers[i].header.data_header.uLength - this->format_containers[i].header.data_header.cLength;
		stats_format[this->header.format_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->format_containers[i].header.stride_header.uLength - this->format_containers[i].header.stride_header.cLength;
		stats_format[this->header.format_offsets[i].data_header.global_key].cost_compressed   += this->format_containers[i].header.data_header.cLength;
		stats_basic[9].cost_uncompressed += this->format_containers[i].getObjectSize();
		stats_basic[9].cost_uncompressed += (S32)this->format_containers[i].header.data_header.uLength - this->format_containers[i].header.data_header.cLength;
		stats_basic[9].cost_uncompressed += (S32)this->format_containers[i].header.stride_header.uLength - this->format_containers[i].header.stride_header.cLength;

	}

	stats_basic[9].cost_compressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	stream.write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));

	std::cerr << "actual: " << last_pos-start_pos << '\t' << "stored: " << this->header.offset_end_of_block << std::endl;

	return(true);
}

}
}