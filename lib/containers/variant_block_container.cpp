#include "variant_block_container.h"

namespace tachyon{
namespace containers{

VariantBlockContainer::VariantBlockContainer() :
	header_(nullptr),
	gt_exp(nullptr)
{

}

VariantBlockContainer::VariantBlockContainer(const global_header_type& header) :
	header_(&header),
	gt_exp(nullptr)
{

}

VariantBlockContainer::VariantBlockContainer(const self_type& other) :
	block_(other.block_),
	compression_manager(other.compression_manager),
	encryption_manager(other.encryption_manager),
	header_(other.header_),
	gt_exp(nullptr)
{
	if(other.gt_exp != nullptr){
		this->gt_exp = new yon_gt_rcd*[this->header_->GetNumberSamples()];
		for(uint32_t i = 0; i < this->header_->GetNumberSamples(); ++i)
			this->gt_exp[i] = other.gt_exp[i];
	}
}

VariantBlockContainer::VariantBlockContainer(self_type&& other) noexcept :
	block_(std::move(other.block_)),
	compression_manager(std::move(other.compression_manager)),
	encryption_manager(std::move(other.encryption_manager)),
	header_(nullptr),
	gt_exp(nullptr)
{
	std::swap(this->gt_exp,  other.gt_exp);
	std::swap(this->header_, other.header_);
}

VariantBlockContainer& VariantBlockContainer::operator=(const self_type& other){
	delete [] this->gt_exp;
	return *this = VariantBlockContainer(other);
}

VariantBlockContainer& VariantBlockContainer::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	this->block_  = std::move(other.block_);
	this->compression_manager = std::move(other.compression_manager);
	this->encryption_manager = std::move(other.encryption_manager);
	this->header_ = nullptr;
	delete [] this->gt_exp; this->gt_exp  = nullptr;
	std::swap(this->gt_exp, other.gt_exp);
	std::swap(this->header_, other.header_);
	return *this;
}

VariantBlockContainer::~VariantBlockContainer(void)
{
	// Do not delete header pointer as this object
	// never owns that data.
	delete [] this->gt_exp;
}

bool VariantBlockContainer::ReadBlock(std::ifstream& stream, block_settings_type& settings){
	return(this->block_.read(stream, settings, *this->header_));
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
	objects->n_loaded_format   = this->block_.load_settings->format_id_global_loaded.size();
	objects->format_containers = new format_interface_type*[this->block_.footer.n_format_streams];

	// Handle Vcf:FORMAT fields.
	if(objects->n_loaded_format){
		for(uint32_t i = 0; i < objects->n_loaded_format; ++i){
			//const uint32_t global_key = this->block_.footer.format_offsets[i].getGlobalKey();
			const uint32_t global_key = this->block_.load_settings->format_id_global_loaded[i];
			const uint32_t local_key  = this->block_.load_settings->format_id_local_loaded[i];
			objects->format_id_loaded.push_back(local_key);
			objects->format_id_loaded_global.push_back(global_key);

			// Evaluate the set-membership of a given global key in the available Format patterns
			// described in the data container footer.
			std::vector<bool> matches = this->block_.FormatPatternSetMembership(global_key);

			if(this->header_->format_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects->format_containers[local_key] = new containers::FormatContainer<int32_t>(this->GetBlock().format_containers[local_key],
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
	objects->n_loaded_info   = this->block_.load_settings->info_id_global_loaded.size();
	objects->info_containers = new info_interface_type*[this->block_.footer.n_info_streams];

	// Handle Vcf:INFO fields.
	if(objects->n_loaded_info){
		for(uint32_t i = 0; i < objects->n_loaded_info; ++i){
			//const uint32_t global_key = this->block_.footer.info_offsets[i].getGlobalKey();
			const uint32_t global_key = this->block_.load_settings->info_id_global_loaded[i];
			const uint32_t local_key  = this->block_.load_settings->info_id_local_loaded[i];
			objects->info_id_loaded.push_back(local_key);
			objects->info_id_loaded_global.push_back(global_key);

			// Evaluate the set-membership of a given global key in the available Info patterns
			// described in the data container footer.
			std::vector<bool> matches = this->block_.InfoPatternSetMembership(global_key);

			if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_INTEGER){
				objects->info_containers[local_key] = new containers::InfoContainer<int32_t>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_STRING ||
					  this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_CHARACTER)
			{
				objects->info_containers[local_key] = new containers::InfoContainer<std::string>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else if(this->header_->info_fields_[global_key].yon_type == YON_VCF_HEADER_FLOAT){
				objects->info_containers[local_key] = new containers::InfoContainer<float>(this->GetBlock().info_containers[local_key], *objects->meta_container, matches);
			} else {
				objects->info_containers[local_key] = new containers::InfoContainer<uint32_t>();
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
	for(uint32_t i = 0; i < objects.meta_container->size(); ++i){
		records[i].is_dirty       = false;
		records[i].is_loaded_meta = true; // todo
		records[i].is_loaded_gt   = objects.loaded_genotypes;
		records[i].id_block       = i;
		records[i].meta = &objects.meta_container->at(i);

		// Check if genotype data has been loaded then checks
		// if this variant site has any genotypes set.
		if(records[i].is_loaded_gt && records[i].meta->HasGT()){
			records[i].gt_i = &objects.genotype_container->at(i);

			if(this->block_.header.controller.has_gt_permuted)
				records[i].gt = records[i].gt_i->GetObjects(*this->block_.gt_ppa);
			else
				records[i].gt = records[i].gt_i->GetObjects(this->header_->GetNumberSamples());

			records[i].gt->Evaluate();
			//records[i].gt->Expand();
		}

		if(objects.n_loaded_info && objects.meta_container->at(i).info_pattern_id >= 0){
			records[i].n_info   = objects.n_loaded_info;
			records[i].info_ids = &this->block_.load_settings->info_patterns_local[objects.meta_container->at(i).info_pattern_id];
			records[i].info     = new PrimitiveContainerInterface*[this->block_.footer.n_info_streams];
			records[i].info_containers = new InfoContainerInterface*[this->block_.footer.n_info_streams];

			// Populate Info data.
			assert(records[i].info_ids != nullptr);
			for(uint32_t j = 0; j < records[i].info_ids->size(); ++j){
				InfoContainerInterface* info_cnt = objects.info_container_map[this->header_->info_fields_[records[i].info_ids->at(j)].id];
				records[i].info_containers[j] = info_cnt;
				records[i].info_hdr.push_back(&this->header_->info_fields_[records[i].info_ids->at(j)]);

				//std::cerr << "INFO " << j << "/" << records[i].info_ids->size() << ": target: " << this->header_->info_fields_[records[i].info_ids->at(j)].id << " container entries = " << info_cnt->size() << std::endl;

				switch(this->header_->info_fields_[records[i].info_ids->at(j)].yon_type){
				case(YON_VCF_HEADER_INTEGER):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<int32_t>*>(info_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_FLOAT):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<float>*>(info_cnt)->at(i);
					break;
				case(YON_VCF_HEADER_STRING):
				case(YON_VCF_HEADER_CHARACTER):
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<std::string>*>(info_cnt)->at(i);
					break;
				default:
					records[i].info[j] = &reinterpret_cast<containers::InfoContainer<uint32_t>*>(info_cnt)->at(0);
					break;
				}
			}
		}

		if(objects.n_loaded_format && objects.meta_container->at(i).format_pattern_id >= 0){
			records[i].n_format   = objects.n_loaded_format;
			records[i].format_ids = &this->block_.load_settings->format_patterns_local[objects.meta_container->at(i).format_pattern_id];
			records[i].fmt        = new PrimitiveGroupContainerInterface*[this->block_.footer.n_format_streams];
			records[i].format_containers = new FormatContainerInterface*[this->block_.footer.n_format_streams];

			// Populate Format data.
			assert(records[i].format_ids != nullptr);
			for(uint32_t j = 0; j < records[i].format_ids->size(); ++j){

				FormatContainerInterface* fmt_cnt = objects.format_container_map[this->header_->format_fields_[records[i].format_ids->at(j)].id];
				records[i].format_containers[j] = fmt_cnt;
				records[i].format_hdr.push_back(&this->header_->format_fields_[records[i].format_ids->at(j)]);

				//std::cerr << "FORMAT " << j << "/" << records[i].format_ids->size() << ": target: " << records[i].format_hdr[j]->id << " container entries = " << fmt_cnt->size() << std::endl;

				switch(this->header_->format_fields_[records[i].format_ids->at(j)].yon_type){
				case(YON_VCF_HEADER_INTEGER):
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<int32_t>*>(fmt_cnt)->at(i);
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
					records[i].fmt[j] = &reinterpret_cast<containers::FormatContainer<uint32_t>*>(fmt_cnt)->at(0);
					break;
				}
			}
		}

		// Todo: These should interpreted during loading as Info and Format
		// fields are.
		if(this->block_.footer.n_filter_streams){
			records[i].n_filter = this->block_.footer.n_filter_streams;
			records[i].filter_ids = &this->block_.footer.filter_patterns[objects.meta_container->at(i).filter_pattern_id].pattern;

			// Populate Filter.
			for(uint32_t j = 0; j < records[i].filter_ids->size(); ++j){
				records[i].filter_hdr.push_back(&this->header_->filter_fields_[records[i].filter_ids->at(j)]);
			}
		}
	}
	return(records);
}

}
}
