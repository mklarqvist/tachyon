#include "variant_block_container.h"

namespace tachyon{
namespace containers{

bool VariantBlockContainer::ReadBlock(std::ifstream& stream, block_settings_type& settings){
	// Todo: slice out target of interest
	// Todo: fix these
	settings.display_static = std::numeric_limits<U32>::max();
	//settings.display_static = (1UL << YON_BLK_GT_INT8) - 1;
	settings.load_static = std::numeric_limits<U32>::max();
	//settings.load_static = (1UL << YON_BLK_GT_INT8) - 1;
	//settings.load_static |= YON_BLK_BV_INFO;

	// Temp: Disable all INFO loading for now.
	settings.load_static &= ~(YON_BLK_BV_INFO);

	if(settings.load_static & YON_BLK_BV_INFO){
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->info_id_loaded.push_back(i);
		}
	} else {
		std::cerr << "ids: " << settings.info_ID_list.size() << std::endl;
		for(U32 i = 0; i < settings.info_ID_list.size(); ++i){
			const int local = this->block_.GetInfoPosition(settings.info_ID_list[i]);
			if(local >= 0) this->info_id_loaded.push_back(local);
			std::cerr << local << std::endl;
		}
	}

	if(settings.load_static & YON_BLK_BV_PPA){
		if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offsets[YON_BLK_PPA].data_header.offset);
			this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_PPA], this->block_.base_containers[YON_BLK_PPA]);
		}
	}

	for(U32 i = YON_BLK_CONTIG; i < YON_BLK_GT_INT8; ++i){
		if(settings.load_static & (1 << i)){
			this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[i], this->block_.base_containers[i]);
		}
	}

	if(settings.load_static & YON_BLK_BV_GT){
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_INT8], this->block_.base_containers[YON_BLK_GT_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT16],   this->block_.base_containers[YON_BLK_GT_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT32],   this->block_.base_containers[YON_BLK_GT_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT64],   this->block_.base_containers[YON_BLK_GT_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT8],  this->block_.base_containers[YON_BLK_GT_S_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT16], this->block_.base_containers[YON_BLK_GT_S_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT32], this->block_.base_containers[YON_BLK_GT_S_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT64], this->block_.base_containers[YON_BLK_GT_S_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT8],  this->block_.base_containers[YON_BLK_GT_N_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT16], this->block_.base_containers[YON_BLK_GT_N_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT32], this->block_.base_containers[YON_BLK_GT_N_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT64], this->block_.base_containers[YON_BLK_GT_N_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_SUPPORT], this->block_.base_containers[YON_BLK_GT_SUPPORT]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_PLOIDY],  this->block_.base_containers[YON_BLK_GT_PLOIDY]);
	}

	// Load all info
	if(this->block_.footer.n_info_streams && (settings.load_static & YON_BLK_BV_INFO)){
		delete [] this->block_.info_containers;
		this->block_.info_containers = new VariantBlock::container_type[this->block_.footer.n_info_streams];
		this->block_.n_info_c_allocated = this->block_.footer.n_info_streams;

		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		//this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
			//std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_info_streams << " -> " << this->block_.footer.info_offsets[i].data_header.global_key << std::endl;
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.info_ID_list.size()){
		std::cerr << "here" << std::endl;
		delete [] this->block_.info_containers;
		this->block_.info_containers = new VariantBlock::container_type[this->block_.footer.n_info_streams];
		this->block_.n_info_c_allocated = this->block_.footer.n_info_streams;

		for(U32 i = 0; i < this->info_id_loaded.size(); ++i){
			this->block_.LoadContainerSeek(stream, this->block_.footer.info_offsets[this->info_id_loaded[i]], this->block_.info_containers[this->info_id_loaded[i]]);
			std::cerr << "Loaded: " << this->info_id_loaded[i] << "/" << this->block_.footer.n_info_streams << " -> " << this->block_.footer.info_offsets[this->info_id_loaded[i]].data_header.global_key << "==" << this->info_id_loaded[i] << std::endl;
		}

	} // end case load_info_ID

	// Load all FORMAT data
	if(this->block_.footer.n_format_streams && (settings.load_static & YON_BLK_BV_FORMAT)){
		delete [] this->block_.format_containers;
		this->block_.format_containers = new VariantBlock::container_type[this->block_.footer.n_format_streams];
		this->block_.n_format_c_allocated = this->block_.footer.n_format_streams;

		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		//this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
			//std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_format_streams << " -> " << this->block_.footer.format_offsets[i].data_header.global_key << std::endl;
		}
		assert(this->block_.end_compressed_data_ == (U64)stream.tellg());
	} // If we have supplied a list of identifiers
	else if(settings.format_ID_list.size()){
		// todo
	} // end case load_info_ID

	stream.seekg(this->block_.end_block_); // seek to end-of-block
	return(true);
}

VariantReaderObjects& VariantBlockContainer::LoadObjects(objects_type& objects, block_settings_type& block_settings){
	// New meta container
	objects.meta_container = new meta_container_type(this->GetBlock());

	// New genotype containers if aplicable
	if(this->HasGenotypes() && (block_settings.load_static & YON_BLK_BV_GT_INT8)){
		objects.genotype_container = new gt_container_type(this->GetBlock(), *objects.meta_container);
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
				objects.format_containers[i] = new containers::FormatContainer<S32>(this->GetBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());

			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects.format_containers[i] = new containers::FormatContainer<std::string>(this->GetBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());
			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects.format_containers[i] = new containers::FormatContainer<float>(this->GetBlock().format_containers[i], *objects.meta_container, matches, this->header_->GetNumberSamples());
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
				objects.info_containers[i] = new containers::InfoContainer<S32>(this->GetBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects.info_containers[i] = new containers::InfoContainer<std::string>(this->GetBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects.info_containers[i] = new containers::InfoContainer<float>(this->GetBlock().info_containers[i], *objects.meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLAG){
				objects.info_containers[i] = new containers::InfoContainer<U32>(true);
			} else {
				objects.info_containers[i] = new containers::InfoContainer<U32>();
			}
			objects.info_field_names.push_back(this->header_->info_fields_[global_key].id);
			objects.info_container_map[this->header_->info_fields_[global_key].id] = objects.info_containers[i];
		}
	}

	return(objects);
}

}
}
