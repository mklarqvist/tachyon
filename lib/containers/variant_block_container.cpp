#include "variant_block_container.h"

namespace tachyon{
namespace containers{

bool VariantBlockContainer::ParseSettings(block_settings_type& settings){
	// Clear previous information (if any).
	this->info_id_global_loaded.clear();
	this->info_id_local_loaded.clear();
	this->info_map_global.clear();
	this->format_id_global_loaded.clear();
	this->format_id_local_loaded.clear();
	this->format_map_global.clear();

	// Todo: if performing genotype annotating then we have to remove
	//       fields that are annotate by tachyon to prevent duplicates
	//       (e.g. AC or AF already existing).
	std::unordered_map<uint32_t, std::string> blocked_list;
	if(settings.annotate_extra){
		for(U32 i = 0; i < YON_GT_ANNOTATE_FIELDS.size(); ++i){
			const YonInfo* info = this->header_->GetInfo(YON_GT_ANNOTATE_FIELDS[i]);
			if(info != nullptr){
				blocked_list[info->idx] = YON_GT_ANNOTATE_FIELDS[i];
				//std::cerr << "Add to blocked list: " << YON_GT_ANNOTATE_FIELDS[i] << "@" << info->idx << std::endl;
			}
		}
	}

	// Parse Info. If all Info containers are loaded then we simply copy
	// the order in which they occur. If we are provided with a vector
	// of target global identifiers we have to first map these to the
	// (possible) local identifiers.
	if(settings.load_static & YON_BLK_BV_INFO){
		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->block_.footer.info_offsets[i].data_header.global_key);
			if(it == blocked_list.end()){
				//std::cerr << "adding not blocked" << std::endl;
				this->info_id_local_loaded.push_back(i);
				this->info_id_global_loaded.push_back(this->block_.footer.info_offsets[i].data_header.global_key);
				this->info_map_global[this->info_id_global_loaded[i]] = i;
			} else {
				//std::cerr << "skipping blocked" << std::endl;
			}
		}
		//std::cerr << this->info_id_local_loaded.size() << "," << this->info_id_global_loaded.size() << std::endl;
	} else {
		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(U32 i = 0; i < settings.info_id_global.size(); ++i){
			// Searches for the global Vcf:INFO idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->block_.footer.info_offsets[i].data_header.global_key);
			if(it == blocked_list.end()){
				const int local = this->block_.GetInfoPosition(settings.info_id_global[i]);
				if(local >= 0){
					local_ids.push_back(local);
					global_ids.push_back(settings.info_id_global[i]);
				}
			}
		}

		if(local_ids.size()){
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
		}
	}

	// Parse Format. If all Format containers are loaded then we simply copy
	// the order in which they occur. If we are provided with a vector
	// of target global identifiers we have to first map these to the
	// (possible) local identifiers.
	if(settings.load_static & YON_BLK_BV_FORMAT){
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->format_id_local_loaded.push_back(i);
			this->format_id_global_loaded.push_back(this->block_.footer.format_offsets[i].data_header.global_key);
			this->format_map_global[this->format_id_global_loaded[i]] = i;
		}
	} else {
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

	return(this->ParseLoadedPatterns(settings));
}

bool VariantBlockContainer::ParseLoadedPatterns(block_settings_type& settings){
	// Clear previous information (if any).
	this->info_patterns_local.clear();
	this->format_patterns_local.clear();
	this->info_patterns_local.resize(this->block_.footer.n_info_patterns);
	this->format_patterns_local.resize(this->block_.footer.n_format_patterns);

	// Iterate over Info patterns.
	if(this->info_id_global_loaded.size()){
		// If all Vcf::INFO fields are desired then return them
		// in the stored order to guarantee bit-exactness. Otherwise
		// return in the order requested.
		if((settings.load_static & YON_BLK_BV_INFO) && settings.annotate_extra == false){
			for(U32 p = 0; p < this->block_.footer.n_info_patterns; ++p){
				this->info_patterns_local[p] = this->block_.footer.info_patterns[p].pattern;
			}
		} else { // Return in requested order.
			for(U32 p = 0; p < this->block_.footer.n_info_patterns; ++p){
				this->info_patterns_local[p] = this->block_.IntersectInfoPatterns(this->info_id_global_loaded, p);
			}
		}
	}

	// Iterate over Format patterns.
	if(this->format_id_global_loaded.size()){
		if(settings.load_static & YON_BLK_BV_FORMAT){
			// If all Vcf::FORMAT fields are desired then return them
			// in the stored order to guarantee bit-exactness. Otherwise
			// return in the order requested.
			for(U32 p = 0; p < this->block_.footer.n_format_patterns; ++p){
				this->format_patterns_local[p] = this->block_.footer.format_patterns[p].pattern;
			}
		} else {
			for(U32 p = 0; p < this->block_.footer.n_format_patterns; ++p){
				this->format_patterns_local[p] = this->block_.IntersectFormatPatterns(this->format_id_global_loaded, p);
			}
		}
	}

	return true;
}

bool VariantBlockContainer::ReadBlock(std::ifstream& stream, block_settings_type& settings){
	// Allocate memory for the Format and Info containers.
	// Info containers.
	delete [] this->block_.info_containers;
	this->block_.info_containers = new VariantBlock::container_type[this->block_.footer.n_info_streams];
	this->block_.n_info_c_allocated = this->block_.footer.n_info_streams;
	// Format containers.
	delete [] this->block_.format_containers;
	this->block_.format_containers = new VariantBlock::container_type[this->block_.footer.n_format_streams];
	this->block_.n_format_c_allocated = this->block_.footer.n_format_streams;

	// Interpret the user-specified block-settings if any. This step converts
	// global index offset values into local offsets and computes new pattern
	// vectors if required. The ordering of the values are according to the
	// input sequence not according to the actual stored order.
	this->ParseSettings(settings);

	// Load the FORMAT:GT (GBPBWT) permutation array.
	if(settings.load_static & YON_BLK_BV_PPA){
		// If there is FORMAT:GT field data available AND that data has
		// been permuted then create a new yon_gt_ppa object to store
		// this data.
		if(this->block_.header.controller.hasGTPermuted && this->block_.header.controller.hasGT){
			stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.offsets[YON_BLK_PPA].data_header.offset);
			this->block_.LoadContainerSeek(stream,
			                               this->block_.footer.offsets[YON_BLK_PPA],
			                               this->block_.base_containers[YON_BLK_PPA]);

			this->block_.gt_ppa = new yon_gt_ppa;
			this->block_.gt_ppa->n_samples = this->header_->GetNumberSamples();
		}
	}

	// Load base meta containers.
	for(U32 i = YON_BLK_CONTIG; i < YON_BLK_GT_INT8; ++i){
		if(settings.load_static & (1 << i)){
			this->block_.LoadContainerSeek(stream,
			                               this->block_.footer.offsets[i],
			                               this->block_.base_containers[i]);
		}
	}

	// Load genotype containers. At the moment, genotype containers cannot be loaded
	// individually by using this wrapper routine. If you wish to load these separately
	// you will have to do so manually.
	if((settings.load_static & YON_BLK_BV_GT) || (settings.load_static & YON_BLK_BV_FORMAT)){
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

	// Load Info containers. Technically there is no difference between the two
	// conditions below in terms of outcome. However, the first case guarantees
	// that data is loaded linearly from disk as this can be guaranteed when loading
	// all available data. There is no such guarntees for the second case.
	if(this->block_.footer.n_info_streams && (settings.load_static & YON_BLK_BV_INFO) && settings.annotate_extra == false){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.info_offsets[0].data_header.offset);

		for(U32 i = 0; i < this->block_.footer.n_info_streams; ++i){
			this->block_.LoadContainer(stream,
			                           this->block_.footer.info_offsets[i],
			                           this->block_.info_containers[i]);
		}
	}
	// If we have a user-supplied list of identifiers parsed above.
	else {
		for(U32 i = 0; i < this->info_id_local_loaded.size(); ++i){
			this->block_.LoadContainerSeek(stream,
			                               this->block_.footer.info_offsets[this->info_id_local_loaded[i]],
			                               this->block_.info_containers[this->info_id_local_loaded[i]]);
		}

	}

	// Load Format containers. Technically there is no difference between the two
	// conditions below in terms of outcome. However, the first case guarantees
	// that data is loaded linearly from disk as this can be guaranteed when loading
	// all available data. There is no such guarntees for the second case.
	if(this->block_.footer.n_format_streams && (settings.load_static & YON_BLK_BV_FORMAT)){
		stream.seekg(this->block_.start_compressed_data_ + this->block_.footer.format_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->block_.footer.n_format_streams; ++i){
			this->block_.LoadContainerSeek(stream, this->block_.footer.format_offsets[i], this->block_.format_containers[i]);
		}
		// At this point the stream should be located at the end-of-block
		// marker as the Format information is stored last.
		assert(this->block_.end_compressed_data_ == (U64)stream.tellg());
	}
	// If we have a user-supplied list of identifiers parsed above.
	else {
		for(U32 i = 0; i < this->format_id_local_loaded.size(); ++i){
			this->block_.LoadContainerSeek(stream, this->block_.footer.format_offsets[this->format_id_local_loaded[i]], this->block_.format_containers[this->format_id_local_loaded[i]]);
		}
	}

	// Seek to end-of-block position.
	stream.seekg(this->block_.end_block_);
	return(true);
}

VariantReaderObjects* VariantBlockContainer::LoadObjects(objects_type* objects, block_settings_type& block_settings){
	// Construct a new high-level meta container using all available loaded
	// core meta data container.
	objects->meta_container = new meta_container_type(this->GetBlock());

	// If the block has genotypes in it and they have been loaded we can construct
	// a high-level genotype container and transform the multifarious internal
	// encodings into a unified framework.
	if(this->HasGenotypes() && (block_settings.load_static & YON_BLK_BV_GT)){
		objects->loaded_genotypes   = true;
		objects->genotype_container = new gt_container_type(this->GetBlock(), *objects->meta_container);
		// Genotype summary object is only allocated when required.
		objects->genotype_summary   = nullptr;
	}

	// Format-specific containers. These have to be allocated as double pointers
	// to avoid memory collisions because they have different intrinsic class members
	// even though they share the same data interface.
	objects->n_loaded_format   = this->format_id_global_loaded.size();
	objects->format_containers = new format_interface_type*[this->block_.footer.n_format_streams];

	// Handle Vcf:FORMAT fields.
	if(objects->n_loaded_format){
		for(U32 i = 0; i < objects->n_loaded_format; ++i){
			//const U32 global_key = this->block_.footer.format_offsets[i].getGlobalKey();
			const U32 global_key = this->format_id_global_loaded[i];
			const U32 local_key  = this->format_id_local_loaded[i];
			objects->format_id_loaded.push_back(local_key);

			// Evaluate the set-membership of a given global key in the available Format patterns
			// described in the data container footer.
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

	// Info-specific containers. These have to be allocated as double pointers
	// to avoid memory collisions because they have different intrinsic class members
	// even though they share the same data interface.
	objects->n_loaded_info   = this->info_id_global_loaded.size();
	objects->info_containers = new info_interface_type*[this->block_.footer.n_info_streams];

	// Handle Vcf:INFO fields.
	if(objects->n_loaded_info){
		for(U32 i = 0; i < objects->n_loaded_info; ++i){
			//const U32 global_key = this->block_.footer.info_offsets[i].getGlobalKey();
			const U32 global_key = this->info_id_global_loaded[i];
			const U32 local_key  = this->info_id_local_loaded[i];
			objects->info_id_loaded.push_back(local_key);

			// Evaluate the set-membership of a given global key in the available Info patterns
			// described in the data container footer.
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
			//records[i].gt->Expand();
		}

		if(objects.n_loaded_info){
			records[i].n_info   = objects.n_loaded_info;
			records[i].info_ids = &this->info_patterns_local[objects.meta_container->at(i).info_pattern_id];
			records[i].info     = new PrimitiveContainerInterface*[this->block_.footer.n_info_streams];
			records[i].info_containers = new InfoContainerInterface*[this->block_.footer.n_info_streams];

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

				//std::cerr << "FORMAT " << j << "/" << records[i].format_ids->size() << ": target: " << records[i].format_hdr[j]->id << " container entries = " << fmt_cnt->size() << std::endl;

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
			}
		}

		// Todo: These should interpreted during loading as Info and Format
		// fields are.
		if(this->block_.footer.n_filter_streams){
			records[i].n_filter = this->block_.footer.n_filter_streams;

			// Populate Filter.
			records[i].filter_ids = &this->block_.footer.filter_patterns[objects.meta_container->at(i).filter_pattern_id].pattern;
			for(U32 j = 0; j < records[i].filter_ids->size(); ++j){
				records[i].filter_hdr.push_back(&this->header_->filter_fields_[records[i].filter_ids->at(j)]);
			}
		}

		// Pointer to this variant block container.
		records[i].parent_container = this;
	}
	return(records);
}

}
}
