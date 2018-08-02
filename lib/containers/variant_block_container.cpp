#include "variant_block_container.h"

namespace tachyon{
namespace containers{

bool VariantBlockContainer::readBlock(std::ifstream& stream, block_settings_type& settings){
	// Todo: slice out target of interest

	// Todo: fix these
	settings.display_static = std::numeric_limits<U32>::max();
	settings.load_static = std::numeric_limits<U32>::max();

	if(settings.display_static & YON_BLK_BV_PPA){
		if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
			this->block_.ppa_manager.header = this->block_.footer.offsets[YON_BLK_PPA];
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offsets[YON_BLK_PPA].data_header.offset);
			stream >> this->block_.ppa_manager;
		}
	}

	for(U32 i = YON_BLK_BV_CONTIG; i < YON_BLK_GT_INT8; ++i){
		if(settings.display_static & (1 << i)){
			this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[i], this->block_.base_containers[i]);
		}
	}

	if(settings.display_static & YON_BLK_BV_GT_INT8){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_INT8], this->block_.base_containers[YON_BLK_GT_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT16], this->block_.base_containers[YON_BLK_GT_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT32], this->block_.base_containers[YON_BLK_GT_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT64], this->block_.base_containers[YON_BLK_GT_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT8], this->block_.base_containers[YON_BLK_GT_S_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT16], this->block_.base_containers[YON_BLK_GT_S_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT32], this->block_.base_containers[YON_BLK_GT_S_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT64], this->block_.base_containers[YON_BLK_GT_S_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_SUPPORT], this->block_.base_containers[YON_BLK_GT_SUPPORT]);
	}
	// Load all info
	if(settings.display_static & YON_BLK_BV_INFO){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		//this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
			std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_info_streams << " -> " << this->block_.footer.info_offsets[i].data_header.global_key << std::endl;
			//this->mapper_.info_container_loaded_.at(i)(i, i, this->block_.footer.info_offsets[i].data_header.global_key, &this->block_.footer.info_offsets[i]);
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.info_ID_list.size()){
		// todo
	} // end case load_info_ID

	// Load all FORMAT data
	if(settings.display_static & YON_BLK_BV_FORMAT){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		//this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
			std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_format_streams << " -> " << this->block_.footer.format_offsets[i].data_header.global_key << std::endl;
			//this->mapper_.format_container_loaded_.at(i)(i, i, this->block_.footer.format_offsets[i].data_header.global_key, &this->block_.footer.format_offsets[i]);
		}
		assert(this->block_.end_compressed_data_ == (U64)stream.tellg());
	} // If we have supplied a list of identifiers
	else if(settings.format_ID_list.size()){
		// todo
	} // end case load_info_ID

	stream.seekg(this->block_.end_block_); // seek to end-of-block
	return(true);
}

VariantReaderObjects& VariantBlockContainer::loadObjects(objects_type& objects, block_settings_type& block_settings){
	// New meta container
	objects.meta_container = new meta_container_type(this->getBlock());

	// New genotype containers if aplicable
	if(this->getBlock().header.controller.hasGT && (block_settings.display_static & YON_BLK_BV_GT_INT8)){
		objects.genotype_container = new gt_container_type(this->getBlock(), *objects.meta_container);
		objects.genotype_summary   = new objects_type::genotype_summary_type(10);
	}

	// FORMAT-specific containers
	// Store as double pointers to avoid memory collisions because
	// FORMAT containers have different intrinsic class members
	objects.n_loaded_format   = this->block_.footer.n_format_streams;
	objects.format_containers = new format_interface_type*[objects.n_loaded_format];

	if(objects.n_loaded_format){
		for(U32 i = 0; i < objects.n_loaded_format; ++i){
			const U32 global_key = this->block_.footer.format_offsets[i].getGlobalKey();

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->block_.FormatPatternSetMembership(global_key);

			if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects.format_containers[i] = new containers::FormatContainer<S32>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());

			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects.format_containers[i] = new containers::FormatContainer<std::string>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());
			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects.format_containers[i] = new containers::FormatContainer<float>(this->getBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());
			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_FLAG){
				std::cerr << utility::timestamp("ERROR") << "Format fields cannot have FLAG fields." << std::endl;
				exit(1);
			} else {
				objects.format_containers[i] = new containers::FormatContainer<U32>;
			}
			objects.format_field_names.push_back(this->header_->format_fields_[global_key].id);
			objects.format_container_map[this->header_->format_fields_[global_key].id] = objects.format_containers[i];
		}
	}

	// INFO-specific containers
	// Store as double pointers to avoid memory collisions because
	// INFO containers have different class members
	objects.n_loaded_info   = this->block_.footer.n_info_streams;
	objects.info_containers = new info_interface_type*[objects.n_loaded_info];

	if(objects.n_loaded_info){
		for(U32 i = 0; i < objects.n_loaded_info; ++i){
			const U32 global_key = this->block_.footer.info_offsets[i].getGlobalKey();

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->block_.InfoPatternSetMembership(global_key);

			if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects.info_containers[i] = new containers::InfoContainer<S32>(this->getBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects.info_containers[i] = new containers::InfoContainer<std::string>(this->getBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects.info_containers[i] = new containers::InfoContainer<float>(this->getBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLAG){
				objects.info_containers[i] = new containers::InfoContainer<U32>(true);
			} else {
				objects.info_containers[i] = new containers::InfoContainer<U32>();
			}
			objects.info_field_names.push_back(this->header_->info_fields_[global_key].id);
			objects.info_container_map[this->header_->info_fields_[global_key].id] = objects.info_containers[i];
		}
	}
	std::cerr << "here after standard" << std::endl;

	// Lazy evaluation of interleaved data.
	yon1_t* records = new yon1_t[objects.meta_container->size()];
	for(U32 i = 0; i < objects.meta_container->size(); ++i){
		records[i].is_dirty = false;
		records[i].is_loaded_meta = true;
		records[i].is_loaded_gt = true;
		records[i].id_block = i;
		records[i].info_ids = &this->block_.footer.info_patterns[objects.meta_container->at(i).info_pattern_id].pattern;
		records[i].info = new PrimitiveContainerInterface*[objects.n_loaded_info];
		records[i].info_containers = new InfoContainerInterface*[objects.n_loaded_info];

		// Populate INFO data.
		for(U32 j = 0; j < records[i].info_ids->size(); ++j){
			InfoContainerInterface* info_cnt = objects.info_container_map[this->header_->info_fields_[records[i].info_ids->at(j)].id];
			records[i].info_containers[j] = info_cnt;
			records[i].info_types.push_back(this->header_->info_fields_[records[i].info_ids->at(j)].yon_type);

			std::cerr << "INFO " << j << "/" << records[i].info_ids->size() << ": target: " << this->header_->info_fields_[records[i].info_ids->at(j)].id << " container entries = " << info_cnt->size() << std::endl;

			switch(this->header_->info_fields_[records[i].info_ids->at(j)].yon_type){
			case(YON_VCF_HEADER_INTEGER):
				records[i].info[j] = &reinterpret_cast<containers::InfoContainer<S32>*>(info_cnt)->at(i);
				break;
			case(YON_VCF_HEADER_FLOAT):
				records[i].info[j] = &reinterpret_cast<containers::InfoContainer<float>*>(info_cnt)->at(i);
				break;
			case(YON_VCF_HEADER_STRING):
			case(YON_VCF_HEADER_CHARACTER):
				records[i].info[j] = &reinterpret_cast<containers::InfoContainer<std::string>*>(info_cnt)->at(i);
				break;
			default:
				records[i].info[j] = &reinterpret_cast<containers::InfoContainer<U32>*>(info_cnt)->at(0);
				break;
			}
			std::cerr << "target entries :" << records[i].info[j]->size() << std::endl;
		}

		std::cerr << "start format" << std::endl;
		records[i].format_ids = &this->block_.footer.format_patterns[objects.meta_container->at(i).format_pattern_id].pattern;
		records[i].fmt = new PrimitiveGroupContainerInterface*[objects.n_loaded_info];
		records[i].format_containers = new FormatContainerInterface*[objects.n_loaded_format];

		std::cerr << "after base format" << std::endl;

		// Populate FORMAT data.
		for(U32 j = 0; j < records[i].format_ids->size(); ++j){
			FormatContainerInterface* fmt_cnt = objects.format_container_map[this->header_->format_fields_[records[i].format_ids->at(j)].id];
			records[i].format_containers[j] = fmt_cnt;
			records[i].format_types.push_back(this->header_->format_fields_[records[i].format_ids->at(j)].yon_type);

			std::cerr << "FORMAT " << j << "/" << records[i].format_ids->size() << ": target: " << this->header_->format_fields_[records[i].format_ids->at(j)].id << " container entries = " << fmt_cnt->size() << std::endl;

			// Todo: if format container is completely empty then skip
			if(fmt_cnt->size() == 0){
				records[i].fmt[j] = nullptr;
				continue;
			}

			switch(this->header_->format_fields_[records[i].format_ids->at(j)].yon_type){
			case(YON_VCF_HEADER_INTEGER):
				records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<S32>*>(fmt_cnt)->at(i);
				break;
			case(YON_VCF_HEADER_FLOAT):
				records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<float>*>(fmt_cnt)->at(i);
				break;
			case(YON_VCF_HEADER_STRING):
			case(YON_VCF_HEADER_CHARACTER):
				records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<std::string>*>(fmt_cnt)->at(i);
				break;
			default:
				std::cerr << "fmt at default" << std::endl;
				records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<U32>*>(fmt_cnt)->at(0);
				break;
			}
			std::cerr << "target entries :" << records[i].fmt[j]->size() << std::endl;
		}

		records[i].parent_container = this;
	}

	std::cerr << "done" << std::endl;

	return(objects);
}

}
}
