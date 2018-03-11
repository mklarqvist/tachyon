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
	for(U32 i = 0; i < this->footer.n_info_streams; ++i)
		this->info_containers[i].reset();

	for(U32 i = 0; i < this->footer.n_format_streams; ++i)
		this->format_containers[i].reset();

	this->header.reset();
	this->footer.reset();
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
	stream >> this->header;
	const U64 start_offset = (U64)stream.tellg();
	stream.seekg(stream.tellg() + this->header.l_offset_footer);
	stream >> this->footer;
	U32 l_return = 0;
	stream.read((char*)&l_return, sizeof(U32));
	U64 eof_marker;
	stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
	assert(eof_marker == constants::TACHYON_BLOCK_EOF);
	const U64 end_of_block = stream.tellg();

	if(settings.importPPA){
		if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
			stream.seekg(start_offset + this->footer.offset_ppa.data_header.offset);
			stream >> this->ppa_manager;
		}
	}

	if(settings.importMetaHot){
		stream.seekg(start_offset + this->footer.offset_hot_meta.data_header.offset);
		this->meta_hot_container.header = this->footer.offset_hot_meta;
		stream >> this->meta_hot_container;
	}

	if(settings.importMetaCold){
		stream.seekg(start_offset + this->footer.offset_cold_meta.data_header.offset);
		this->meta_cold_container.header = this->footer.offset_cold_meta;
		stream >> this->meta_cold_container;
	}

	if(settings.importGT && this->footer.offset_gt_rle.data_header.offset != -1){
		stream.seekg(start_offset + this->footer.offset_gt_rle.data_header.offset);
		this->gt_rle_container.header = this->footer.offset_gt_rle;
		stream >> this->gt_rle_container;
	}

	if(settings.importGTSimple && this->footer.offset_gt_simple.data_header.offset != -1){
		stream.seekg(start_offset + this->footer.offset_gt_simple.data_header.offset);
		this->gt_simple_container.header = this->footer.offset_gt_simple;
		stream >> this->gt_simple_container;
	}

	if((settings.importGTSimple || settings.importGT) && this->footer.offset_gt_helper.data_header.offset != -1){
		stream.seekg(start_offset + this->footer.offset_gt_helper.data_header.offset);
		this->gt_support_data_container.header = this->footer.offset_gt_helper;
		stream >> this->gt_support_data_container;
	}

	if((settings.importMetaHot || settings.importMetaCold) && this->footer.offset_meta_info_id.data_header.offset != -1){
		stream.seekg(start_offset + this->footer.offset_meta_info_id.data_header.offset);
		this->meta_info_map_ids.header = this->footer.offset_meta_info_id;
		stream >> this->meta_info_map_ids;
		this->meta_filter_map_ids.header = this->footer.offset_meta_filter_id;
		stream >> this->meta_filter_map_ids;
		this->meta_format_map_ids.header = this->footer.offset_meta_format_id;
		stream >> this->meta_format_map_ids;
	}

	// Load all info
	if(settings.importInfoAll && this->footer.n_info_streams > 0){
		stream.seekg(start_offset + this->footer.info_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_info_streams; ++i){
			this->info_containers[i].header = this->footer.info_offsets[i];
			stream >> this->info_containers[i];
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.info_ID_list.size() > 0) {
		core::SettingsMap map_entry;
		// Cycle over all the keys we are interested in
		U32 iterator_index = 0;
		for(U32 i = 0; i < settings.info_ID_list.size(); ++i){
			// Cycle over all available streams in this block
			for(U32 j = 0; j < this->footer.n_info_streams; ++j){
				// If there is a match
				// Push back field into map
				if(this->footer.info_offsets[j].data_header.global_key == settings.info_ID_list[i]){
					settings.load_info_ID_loaded.push_back(
							core::SettingsMap(
									iterator_index++,                   // iterator value
									j,                                  // local index id
									&this->footer.info_offsets[j]) // offset
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
			this->info_containers[settings.load_info_ID_loaded[i].iterator_index].header = this->footer.info_offsets[settings.load_info_ID_loaded[i].iterator_index];
			stream >> this->info_containers[settings.load_info_ID_loaded[i].iterator_index];
		}
	} // end case load_info_ID

	if(settings.importFormatAll && this->footer.n_format_streams){
		stream.seekg(start_offset + this->footer.format_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_format_streams; ++i){
			this->format_containers[i].header = this->footer.format_offsets[i];
			stream >> this->format_containers[i];
			//std::cerr << "loaded: " << this->index_entry.format_offsets[i].global_key << '\t' << this->format_containers[i].header.cLength << std::endl;
		}
	}

	stream.seekg(end_of_block);
	return(true);
}

const U64 VariantBlock::__determineCompressedSize(void) const{
	U64 total = 0;
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		total += this->ppa_manager.getObjectSize();

	total += this->meta_hot_container.getObjectSize();
	total += this->meta_cold_container.getObjectSize();
	total += this->gt_rle_container.getObjectSize();
	total += this->gt_simple_container.getObjectSize();
	total += this->gt_support_data_container.getObjectSize();
	total += this->meta_info_map_ids.getObjectSize();
	total += this->meta_filter_map_ids.getObjectSize();
	total += this->meta_format_map_ids.getObjectSize();
	for(U32 i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].getObjectSize();
	for(U32 i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].getObjectSize();

	return(total);
}

bool VariantBlock::write(std::ofstream& stream,
                         import_stats_type& stats_basic,
                         import_stats_type& stats_info,
                         import_stats_type& stats_format)
{
	U64 last_pos = stream.tellp();
	this->header.l_offset_footer = this->__determineCompressedSize();
	stream << this->header;
	const U64 start_pos = stream.tellp();
	stats_basic[0].cost_uncompressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
		this->footer.offset_ppa.data_header.offset = stream.tellp() - start_pos;
		stream << this->ppa_manager;
		stats_basic[1].cost_compressed   += (U64)stream.tellp() - last_pos;
		stats_basic[1].cost_uncompressed += this->ppa_manager.getObjectSize();
		last_pos = stream.tellp();
	}

	this->__updateHeader(this->footer.offset_hot_meta, this->meta_hot_container, stream.tellp() - start_pos);
	stream << this->meta_hot_container;
	stats_basic[2].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[2].cost_uncompressed += this->meta_hot_container.getObjectSize();
	stats_basic[2].cost_uncompressed += (S32)this->meta_hot_container.header.data_header.uLength - this->meta_hot_container.header.data_header.cLength;
	stats_basic[2].cost_uncompressed += (S32)this->meta_hot_container.header.stride_header.uLength - this->meta_hot_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_cold_meta, this->meta_cold_container, stream.tellp() - start_pos);
	stream << this->meta_cold_container;
	stats_basic[3].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[3].cost_uncompressed += this->meta_cold_container.getObjectSize();
	stats_basic[3].cost_uncompressed += (S32)this->meta_cold_container.header.data_header.uLength - this->meta_cold_container.header.data_header.cLength;
	stats_basic[3].cost_uncompressed += (S32)this->meta_cold_container.header.stride_header.uLength - this->meta_cold_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_gt_rle, this->gt_rle_container, stream.tellp() - start_pos);
	stream << this->gt_rle_container;
	stats_basic[4].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[4].cost_uncompressed += this->gt_rle_container.getObjectSize();
	stats_basic[4].cost_uncompressed += (S32)this->gt_rle_container.header.data_header.uLength - this->gt_rle_container.header.data_header.cLength;
	stats_basic[4].cost_uncompressed += (S32)this->gt_rle_container.header.stride_header.uLength - this->gt_rle_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_gt_simple, this->gt_simple_container, stream.tellp() - start_pos);
	stream << this->gt_simple_container;
	stats_basic[5].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[5].cost_uncompressed += this->gt_simple_container.getObjectSize();
	stats_basic[5].cost_uncompressed += (S32)this->gt_simple_container.header.data_header.uLength - this->gt_simple_container.header.data_header.cLength;
	stats_basic[5].cost_uncompressed += (S32)this->gt_simple_container.header.stride_header.uLength - this->gt_simple_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_gt_helper, this->gt_support_data_container, stream.tellp() - start_pos);
	stream << this->gt_support_data_container;
	stats_basic[6].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[6].cost_uncompressed += this->gt_support_data_container.getObjectSize();
	stats_basic[6].cost_uncompressed += (S32)this->gt_support_data_container.header.data_header.uLength - this->gt_support_data_container.header.data_header.cLength;
	stats_basic[6].cost_uncompressed += (S32)this->gt_support_data_container.header.stride_header.uLength - this->gt_support_data_container.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_meta_info_id, this->meta_info_map_ids, stream.tellp() - start_pos);
	stream << this->meta_info_map_ids;
	stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[7].cost_uncompressed += this->meta_info_map_ids.getObjectSize();
	stats_basic[7].cost_uncompressed += (S32)this->meta_info_map_ids.header.data_header.uLength - this->meta_info_map_ids.header.data_header.cLength;
	stats_basic[7].cost_uncompressed += (S32)this->meta_info_map_ids.header.stride_header.uLength - this->meta_info_map_ids.header.stride_header.cLength;
	last_pos = stream.tellp();



	this->__updateHeader(this->footer.offset_meta_filter_id, this->meta_filter_map_ids, stream.tellp() - start_pos);
	stream << this->meta_filter_map_ids;
	stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[7].cost_uncompressed += this->meta_filter_map_ids.getObjectSize();
	stats_basic[7].cost_uncompressed += (S32)this->meta_filter_map_ids.header.data_header.uLength - this->meta_filter_map_ids.header.data_header.cLength;
	stats_basic[7].cost_uncompressed += (S32)this->meta_filter_map_ids.header.stride_header.uLength - this->meta_filter_map_ids.header.stride_header.cLength;
	last_pos = stream.tellp();

	this->__updateHeader(this->footer.offset_meta_format_id, this->meta_format_map_ids, stream.tellp() - start_pos);
	stream << this->meta_format_map_ids;
	stats_basic[7].cost_compressed   += (U64)stream.tellp() - last_pos;
	stats_basic[7].cost_uncompressed += this->meta_format_map_ids.getObjectSize();
	stats_basic[7].cost_uncompressed += (S32)this->meta_format_map_ids.header.data_header.uLength - this->meta_format_map_ids.header.data_header.cLength;
	stats_basic[7].cost_uncompressed += (S32)this->meta_format_map_ids.header.stride_header.uLength - this->meta_format_map_ids.header.stride_header.cLength;
	last_pos = stream.tellp();

	for(U32 i = 0; i < this->footer.n_info_streams; ++i){
		this->__updateHeader(this->footer.info_offsets[i], this->info_containers[i], stream.tellp() - start_pos);
		stream << this->info_containers[i];
		stats_info[this->footer.info_offsets[i].data_header.global_key].cost_uncompressed += this->info_containers[i].header.data_header.uLength;
		stats_info[this->footer.info_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->info_containers[i].header.data_header.uLength - this->info_containers[i].header.data_header.cLength;
		stats_info[this->footer.info_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->info_containers[i].header.stride_header.uLength - this->info_containers[i].header.stride_header.cLength;
		stats_info[this->footer.info_offsets[i].data_header.global_key].cost_compressed   += this->info_containers[i].header.data_header.cLength;
		stats_basic[8].cost_uncompressed += this->info_containers[i].getObjectSize();
		stats_basic[8].cost_uncompressed += (S32)this->info_containers[i].header.data_header.uLength - this->info_containers[i].header.data_header.cLength;
		stats_basic[8].cost_uncompressed += (S32)this->info_containers[i].header.stride_header.uLength - this->info_containers[i].header.stride_header.cLength;
	}

	stats_basic[8].cost_compressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	for(U32 i = 0; i < this->footer.n_format_streams; ++i){
		this->__updateHeader(this->footer.format_offsets[i], this->format_containers[i], stream.tellp() - start_pos);
		stream << this->format_containers[i];
		stats_format[this->footer.format_offsets[i].data_header.global_key].cost_uncompressed += this->format_containers[i].header.data_header.uLength;
		stats_format[this->footer.format_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->format_containers[i].header.data_header.uLength - this->format_containers[i].header.data_header.cLength;
		stats_format[this->footer.format_offsets[i].data_header.global_key].cost_uncompressed += (S32)this->format_containers[i].header.stride_header.uLength - this->format_containers[i].header.stride_header.cLength;
		stats_format[this->footer.format_offsets[i].data_header.global_key].cost_compressed   += this->format_containers[i].header.data_header.cLength;
		stats_basic[9].cost_uncompressed += this->format_containers[i].getObjectSize();
		stats_basic[9].cost_uncompressed += (S32)this->format_containers[i].header.data_header.uLength - this->format_containers[i].header.data_header.cLength;
		stats_basic[9].cost_uncompressed += (S32)this->format_containers[i].header.stride_header.uLength - this->format_containers[i].header.stride_header.cLength;
	}

	stats_basic[9].cost_compressed += (U64)stream.tellp() - last_pos;
	last_pos = stream.tellp();

	// writing footer
	assert(this->header.l_offset_footer == stream.tellp() - start_pos);
	stream << this->footer;
	const U32 return_offset = stream.tellp() - last_pos;
	// Write negative offset to start of offsets
	stream.write(reinterpret_cast<const char*>(&return_offset), sizeof(U32));
	// Write EOB
	stream.write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));

	return(true);
}

}
}
