#include "variant_block_container.h"

namespace tachyon{
namespace containers{

bool VariantBlockContainer::ParseSettings(block_settings_type& settings){
	this->info_id_global_loaded.clear();
	this->info_id_local_loaded.clear();
	this->info_map_global.clear();
	this->format_id_global_loaded.clear();
	this->format_id_local_loaded.clear();
	this->format_map_global.clear();

	// Parse Info.
	if(settings.load_static & YON_BLK_BV_INFO){
		//std::cerr << "add all info" << std::endl;
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->info_id_local_loaded.push_back(i);
			this->info_id_global_loaded.push_back(this->block_.footer.info_offsets[i].data_header.global_key);
			this->info_map_global[this->info_id_global_loaded[i]] = i;
		}
	} else {
		//std::cerr << "info ids: " << settings.info_id_global.size() << std::endl;

		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(U32 i = 0; i < settings.info_id_global.size(); ++i){
			// Searches for the global Vcf:INFO idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			const int local = this->block_.GetInfoPosition(settings.info_id_global[i]);
			if(local >= 0){
				local_ids.push_back(local);
				global_ids.push_back(settings.info_id_global[i]);
			}
			//std::cerr << "info local: " << local << std::endl;
		}

		if(local_ids.size()){
			//std::cerr << "deduping local info" << std::endl;
			// Dedupe vectors. This prevents multiple parsings of the same
			// target data container as this is illegal.
			for(U32 i = 0; i < local_ids.size(); ++i){
				map_type::const_iterator it = this->info_map_global.find(global_ids[i]);
				if(it == this->info_map_global.end()){
					this->info_id_local_loaded.push_back(local_ids[i]);
					this->info_id_global_loaded.push_back(global_ids[i]);
					this->info_map_global[global_ids[i]] = i;
				}
			}
			//std::cerr << "loaded info: " << this->info_id_local_loaded.size() << std::endl;
		}
	}

	// Parse Format.
	if(settings.load_static & YON_BLK_BV_FORMAT){
		//std::cerr << "add all format: " << this->block_.footer.n_format_streams << std::endl;
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			//std::cerr << "format adding: " << this->block_.footer.format_offsets[i].data_header.global_key << "->" << this->header_->format_fields_[this->block_.footer.format_offsets[i].data_header.global_key].id << std::endl;
			this->format_id_local_loaded.push_back(i);
			this->format_id_global_loaded.push_back(this->block_.footer.format_offsets[i].data_header.global_key);
			this->format_map_global[this->format_id_global_loaded[i]] = i;
		}
	} else {
		//std::cerr << "format ids: " << settings.format_id_global.size() << std::endl;

		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(U32 i = 0; i < settings.format_id_global.size(); ++i){
			// Searches for the global Vcf:FORMAT idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			const int local = this->block_.GetFormatPosition(settings.format_id_global[i]);
			if(local >= 0){
				local_ids.push_back(local);
				global_ids.push_back(settings.format_id_global[i]);
			}
			//std::cerr << "local: " << local << std::endl;
		}

		if(local_ids.size()){
			// Dedupe vectors. This prevents multiple parsings of the same
			// target data container as this is illegal.
			for(U32 i = 0; i < local_ids.size(); ++i){
				map_type::const_iterator it = this->format_map_global.find(global_ids[i]);
				if(it == this->format_map_global.end()){
					this->format_id_local_loaded.push_back(local_ids[i]);
					this->format_id_global_loaded.push_back(global_ids[i]);
					this->format_map_global[global_ids[i]] = i;
				}
			}
		}
	}

	return(this->ParseLoadedPatterns());
}

bool VariantBlockContainer::ParseLoadedPatterns(){
	this->info_patterns_local.clear();
	this->format_patterns_local.clear();
	this->info_patterns_local.resize(this->block_.footer.n_info_patterns);
	this->format_patterns_local.resize(this->block_.footer.n_format_patterns);

	// Iterate over Info patterns.
	if(this->info_id_global_loaded.size()){
		for(U32 p = 0; p < this->block_.footer.n_info_patterns; ++p){
			this->info_patterns_local[p] = this->block_.IntersectInfoPatterns(this->info_id_global_loaded, p);
		}
	}

	// Iterate over Format patterns.
	if(this->format_id_global_loaded.size()){
		for(U32 p = 0; p < this->block_.footer.n_format_patterns; ++p){
			this->format_patterns_local[p] = this->block_.IntersectFormatPatterns(this->format_id_global_loaded, p);
		}
	}

	return true;
}

bool VariantBlockContainer::ReadBlock(std::ifstream& stream, block_settings_type& settings){
	// Todo: slice out target of interest
	// Todo: fix these
	//settings.display_static = std::numeric_limits<U32>::max();
	//settings.display_static = (1UL << YON_BLK_GT_INT8) - 1;
	//settings.load_static = std::numeric_limits<U32>::max();
	//settings.load_static = (1UL << YON_BLK_GT_INT8) - 1;
	//settings.load_static |= YON_BLK_BV_INFO;

	// Temp: Disable all INFO and FORMAT loading.
	//settings.load_static &= ~(YON_BLK_BV_PPA);
	//settings.load_static &= ~(YON_BLK_BV_INFO);
	//settings.load_static &= ~(YON_BLK_BV_FORMAT);
	//settings.load_static &= ~(YON_BLK_BV_GT);
	std::cerr << "Flags: " << std::bitset<32>(settings.load_static) << std::endl;

	// Restore
	delete [] this->block_.info_containers;
	this->block_.info_containers = new VariantBlock::container_type[this->block_.footer.n_info_streams];
	this->block_.n_info_c_allocated = this->block_.footer.n_info_streams;
	delete [] this->block_.format_containers;
	this->block_.format_containers = new VariantBlock::container_type[this->block_.footer.n_format_streams];
	this->block_.n_format_c_allocated = this->block_.footer.n_format_streams;

	// parse settings
	this->ParseSettings(settings);

	if(settings.load_static & YON_BLK_BV_PPA){
		// If there is FORMAT:GT field data available AND that data has
		// been permuted then create a new yon_gt_ppa object to store
		// this data.
		if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offsets[YON_BLK_PPA].data_header.offset);
			this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_PPA], this->block_.base_containers[YON_BLK_PPA]);

			this->block_.gt_ppa = new yon_gt_ppa;
			this->block_.gt_ppa->n_samples = this->header_->GetNumberSamples();
		}
	}

	for(U32 i = YON_BLK_CONTIG; i < YON_BLK_GT_INT8; ++i){
		if(settings.load_static & (1 << i)){
			this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[i], this->block_.base_containers[i]);
		}
	}

	if(settings.load_static & YON_BLK_BV_GT){
		//std::cerr << "loading genotypes" << std::endl;
		this->loaded_genotypes = true;
		this->block_.LoadContainerSeek(stream, this->block_.footer.offsets[YON_BLK_GT_INT8], this->block_.base_containers[YON_BLK_GT_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT16],    this->block_.base_containers[YON_BLK_GT_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT32],    this->block_.base_containers[YON_BLK_GT_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_INT64],    this->block_.base_containers[YON_BLK_GT_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT8],   this->block_.base_containers[YON_BLK_GT_S_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT16],  this->block_.base_containers[YON_BLK_GT_S_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT32],  this->block_.base_containers[YON_BLK_GT_S_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_S_INT64],  this->block_.base_containers[YON_BLK_GT_S_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT8],   this->block_.base_containers[YON_BLK_GT_N_INT8]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT16],  this->block_.base_containers[YON_BLK_GT_N_INT16]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT32],  this->block_.base_containers[YON_BLK_GT_N_INT32]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_N_INT64],  this->block_.base_containers[YON_BLK_GT_N_INT64]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_SUPPORT],  this->block_.base_containers[YON_BLK_GT_SUPPORT]);
		this->block_.LoadContainer(stream, this->block_.footer.offsets[YON_BLK_GT_PLOIDY],   this->block_.base_containers[YON_BLK_GT_PLOIDY]);
	}

	// Load all info
	if(this->block_.footer.n_info_streams && (settings.load_static & YON_BLK_BV_INFO)){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		//this->mapper_.info_container_loaded_.resize(this->block_.footer.n_info_streams);
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.info_offsets[i], this->block_.info_containers[i]);
			//std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_info_streams << " -> " << this->block_.footer.info_offsets[i].data_header.global_key << std::endl;
		}
	}
	// If we have supplied a list of identifiers
	else {
		//std::cerr << "here load selected info: " << this->info_id_local_loaded.size() << std::endl;

		for(U32 i = 0; i < this->info_id_local_loaded.size(); ++i){
			//std::cerr << "Loaded: " << i << "/" << this->info_id_local_loaded.size() << " -> " << this->info_id_local_loaded[i] << "/" << this->block_.footer.n_info_streams << " -> " << this->block_.footer.info_offsets[this->info_id_local_loaded[i]].data_header.global_key << "==" << this->info_id_local_loaded[i] << std::endl;
			//std::cerr << i << "/" << this->block_.n_info_c_allocated << " and " << std::endl;
			this->block_.LoadContainerSeek(stream, this->block_.footer.info_offsets[this->info_id_local_loaded[i]], this->block_.info_containers[this->info_id_local_loaded[i]]);
		}

	} // end case load_info_ID

	// Load all FORMAT data
	//std::cerr << "Match all format: " << std::bitset<32>(settings.load_static & YON_BLK_BV_FORMAT) << std::endl;
	if(this->block_.footer.n_format_streams && (settings.load_static & YON_BLK_BV_FORMAT)){
		//std::cerr << "load all format" << std::endl;

		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		//this->mapper_.format_container_loaded_.resize(this->block_.footer.n_format_streams);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.LoadContainer(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
			//std::cerr << "Loaded: " << i << "/" << this->block_.footer.n_format_streams << " -> " << this->block_.footer.format_offsets[i].data_header.global_key << std::endl;
		}
		assert(this->block_.end_compressed_data_ == (U64)stream.tellg());
	} // If we have supplied a list of identifiers
	else {
		//std::cerr << "here format: " << this->format_id_local_loaded.size() << std::endl;

		for(U32 i = 0; i < this->format_id_local_loaded.size(); ++i){
			this->block_.LoadContainerSeek(stream, this->block_.footer.format_offsets[this->format_id_local_loaded[i]], this->block_.format_containers[this->format_id_local_loaded[i]]);
			//std::cerr << "Loaded: " << this->format_id_local_loaded[i] << "/" << this->block_.footer.n_format_streams << " -> " << this->block_.footer.format_offsets[this->format_id_local_loaded[i]].data_header.global_key << "==" << this->info_id_local_loaded[i] << std::endl;
		}
	} // end case load_info_ID

	stream.seekg(this->block_.end_block_); // seek to end-of-block
	//std::cerr << "done loading" << std::endl;
	return(true);
}

VariantReaderObjects* VariantBlockContainer::LoadObjects(objects_type* objects, block_settings_type& block_settings){
	// New meta container.
	objects->meta_container = new meta_container_type(this->GetBlock());

	// New genotype containers if applicable.
	if(this->HasGenotypes() && (block_settings.load_static & YON_BLK_BV_GT)){
		objects->loaded_genotypes   = true;
		objects->genotype_container = new gt_container_type(this->GetBlock(), *objects->meta_container);
		objects->genotype_summary   = new objects_type::genotype_summary_type(10);
	}

	// FORMAT-specific containers.
	// Store as double pointers to avoid memory collisions because
	// FORMAT containers have different intrinsic class members.
	objects->n_loaded_format   = this->format_id_global_loaded.size();
	objects->format_containers = new format_interface_type*[this->block_.footer.n_format_streams];

	if(objects->n_loaded_format){
		for(U32 i = 0; i < objects->n_loaded_format; ++i){
			//const U32 global_key = this->block_.footer.format_offsets[i].getGlobalKey();
			const U32 global_key = this->format_id_global_loaded[i];
			const U32 local_key  = this->format_id_local_loaded[i];
			objects->format_id_loaded.push_back(local_key);

			//std::cerr << "format: " << i << "/" << objects->n_loaded_format << " -> " << global_key << ","
			//          << local_key << " == " << this->header_->format_fields_[global_key].id << ":"
			//          << this->header_->format_fields_[global_key].type << std::endl;

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns.
			std::vector<bool> matches = this->block_.FormatPatternSetMembership(global_key);

			if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects->format_containers[local_key] = new containers::FormatContainer<S32>(this->GetBlock().format_containers[local_key],
				                                                                             *objects->meta_container,
				                                                                             matches,
				                                                                             this->header_->GetNumberSamples());

			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects->format_containers[local_key] = new containers::FormatContainer<std::string>(this->GetBlock().format_containers[local_key],
				                                                                                     *objects->meta_container,
				                                                                                     matches,
				                                                                                     this->header_->GetNumberSamples());
			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects->format_containers[local_key] = new containers::FormatContainer<float>(this->GetBlock().format_containers[local_key],
				                                                                               *objects->meta_container,
				                                                                               matches,
				                                                                               this->header_->GetNumberSamples());
			} else if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_FLAG){
				std::cerr << utility::timestamp("ERROR") << "Format fields cannot have FLAG fields." << std::endl;
				exit(1);
			}
			objects->format_container_map[this->header_->format_fields_[global_key].id] = objects->format_containers[local_key];
		}
	}



	// INFO-specific containers
	// Store as double pointers to avoid memory collisions because
	// INFO containers have different class members
	objects->n_loaded_info   = this->info_id_global_loaded.size();
	objects->info_containers = new info_interface_type*[this->block_.footer.n_info_streams];

	if(objects->n_loaded_info){
		for(U32 i = 0; i < objects->n_loaded_info; ++i){
			//const U32 global_key = this->block_.footer.info_offsets[i].getGlobalKey();
			const U32 global_key = this->info_id_global_loaded[i];
			const U32 local_key  = this->info_id_local_loaded[i];
			objects->info_id_loaded.push_back(local_key);

			//std::cerr << "building info: " << i << "/" << objects->n_loaded_info << " -> " << global_key << "," << local_key << " == " << this->header_->info_fields_[global_key].id << std::endl;

			// Pattern matches of GLOBAL in LOCAL
			// This evaluated the boolean set-membership of GLOBAL key in the FORMAT patterns
			std::vector<bool> matches = this->block_.InfoPatternSetMembership(global_key);

			if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects->info_containers[local_key] = new containers::InfoContainer<S32>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects->info_containers[local_key] = new containers::InfoContainer<std::string>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects->info_containers[local_key] = new containers::InfoContainer<float>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else {
				objects->info_containers[local_key] = new containers::InfoContainer<U32>();
			}
			objects->info_container_map[this->header_->info_fields_[global_key].id] = objects->info_containers[local_key];
		}
	}

	return(objects);
}

yon1_t* VariantBlockContainer::LazyEvaluate(objects_type& objects){
	// Lazy evaluation of interleaved data.
	yon1_t* records = new yon1_t[objects.meta_container->size()];

	// Iterate over the sites described in the meta container.
	for(U32 i = 0; i < objects.meta_container->size(); ++i){
		records[i].is_dirty       = false;
		records[i].is_loaded_meta = true; // todo
		records[i].is_loaded_gt   = objects.loaded_genotypes;
		records[i].id_block       = i;
		records[i].meta = &objects.meta_container->at(i);

		// Check if genotype data has been loaded then checks
		// if this variant site has any genotypes set.
		if(records[i].is_loaded_gt && records[i].meta->HasGT()){
			records[i].gt_i = &objects.genotype_container->at(i);

			if(this->block_.header.controller.hasGTPermuted)
				records[i].gt = records[i].gt_i->GetObjects(*this->block_.gt_ppa);
			else
				records[i].gt = records[i].gt_i->GetObjects(this->header_->GetNumberSamples());

			records[i].gt->Evaluate();
			records[i].gt->Expand();
		}

		if(objects.n_loaded_info){
			records[i].n_info   = objects.n_loaded_info;
			records[i].info_ids = &this->info_patterns_local[objects.meta_container->at(i).info_pattern_id];
			records[i].info     = new PrimitiveContainerInterface*[this->block_.footer.n_info_streams];
			records[i].info_containers = new InfoContainerInterface*[this->block_.footer.n_info_streams];
			//std::cerr << "here in loaded info: " << records[i].info_ids->size() << std::endl;

			// Populate Info data.
			for(U32 j = 0; j < records[i].info_ids->size(); ++j){
				InfoContainerInterface* info_cnt = objects.info_container_map[this->header_->info_fields_[records[i].info_ids->at(j)].id];
				records[i].info_containers[j] = info_cnt;
				records[i].info_hdr.push_back(&this->header_->info_fields_[records[i].info_ids->at(j)]);

				//std::cerr << "INFO " << j << "/" << records[i].info_ids->size() << ": target: " << this->header_->info_fields_[records[i].info_ids->at(j)].id << " container entries = " << info_cnt->size() << std::endl;

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
				//std::cerr << "target entries :" << records[i].info[j]->size() << std::endl;
			}
		}

		if(objects.n_loaded_format){
			records[i].n_format   = objects.n_loaded_format;
			records[i].format_ids = &this->format_patterns_local[objects.meta_container->at(i).format_pattern_id];
			records[i].fmt        = new PrimitiveGroupContainerInterface*[this->block_.footer.n_format_streams];
			records[i].format_containers = new FormatContainerInterface*[this->block_.footer.n_format_streams];

			// Populate Format data.
			for(U32 j = 0; j < records[i].format_ids->size(); ++j){
				FormatContainerInterface* fmt_cnt = objects.format_container_map[this->header_->format_fields_[records[i].format_ids->at(j)].id];
				records[i].format_containers[j] = fmt_cnt;
				records[i].format_hdr.push_back(&this->header_->format_fields_[records[i].format_ids->at(j)]);

				//std::cerr << "FORMAT " << j << "/" << records[i].format_ids->size() << ": target: " << this->header_->format_fields_[records[i].format_ids->at(j)].id << " container entries = " << fmt_cnt->size() << std::endl;

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
					std::cerr << "fmt at default: illegal primitive in assignment" << std::endl;
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<U32>*>(fmt_cnt)->at(0);
					break;
				}
				//if(this->header_->format_fields_[records[i].format_ids->at(j)].id != "GT")
				//	std::cerr << "target entries :" << records[i].fmt[j]->size() << std::endl;
			}
		}

		if(this->block_.footer.n_filter_streams){ // Todo
			records[i].n_filter = this->block_.footer.n_filter_streams;

			// Populate Filter.
			records[i].filter_ids = &this->block_.footer.filter_patterns[objects.meta_container->at(i).filter_pattern_id].pattern;
			for(U32 j = 0; j < records[i].filter_ids->size(); ++j){
				records[i].filter_hdr.push_back(&this->header_->filter_fields_[records[i].filter_ids->at(j)]);
			}
		}

		records[i].parent_container = this;
	}

	//std::cerr << "done" << std::endl;
	return(records);
}

}
}
