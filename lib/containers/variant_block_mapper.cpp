#include "variant_block_mapper.h"

namespace tachyon{
namespace containers{

VariantBlockMapper::VariantBlockMapper() :
	n_format_fields(0),
	n_info_fields(0)
{
}

VariantBlockMapper::VariantBlockMapper(const size_t n_format_fields, const size_t n_info_fields) :
	n_format_fields(n_format_fields),
	n_info_fields(n_info_fields)
{
}

VariantBlockMapper::~VariantBlockMapper(){}

bool VariantBlockMapper::build(const block_footer_type& block_footer){
	this->format_container_global_.clear();
	this->format_container_global_.resize(this->n_format_fields);
	this->format_container_local_.clear();
	this->format_container_local_.resize(block_footer.n_format_streams);

	for(U32 i = 0; i < block_footer.n_format_streams; ++i){
		// Set global -> local mapping -> loaded mapping
		this->format_container_global_[block_footer.format_offsets[i].data_header.global_key](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
		// Set local -> global mapping -> loaded mapping
		this->format_container_local_[i](i, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
	}

	this->info_container_global_.clear();
	this->info_container_global_.resize(this->n_info_fields);
	this->info_container_local_.clear();
	this->info_container_local_.resize(block_footer.n_info_streams);

	for(U32 i = 0; i < block_footer.n_info_streams; ++i){
		// Set global -> local mapping -> loaded mapping
		this->info_container_global_[block_footer.info_offsets[i].data_header.global_key](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
		// Set local -> global mapping -> loaded mapping
		this->info_container_local_[i](i, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
	}


	return true;
}

bool VariantBlockMapper::build(const std::vector<U32>& info_keys, const std::vector<U32>& format_keys, const block_footer_type& block_footer){
	this->format_container_global_.clear();
	this->format_container_global_.resize(this->n_format_fields);
	this->format_container_local_.clear();
	this->format_container_local_.resize(block_footer.n_format_streams);
	this->format_container_loaded_.clear();

	for(U32 i = 0; i < block_footer.n_format_streams; ++i){
		// Set global -> local mapping -> loaded mapping
		this->format_container_global_[block_footer.format_offsets[i].data_header.global_key](-1, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
		// Set local -> global mapping -> loaded mapping
		this->format_container_local_[i](-1, i, block_footer.format_offsets[i].data_header.global_key, &block_footer.format_offsets[i]);
	}

	for(U32 i = 0; i < format_keys.size(); ++i){
		this->format_container_global_[format_keys[i]].load_order_index = i;
		this->format_container_local_[this->format_container_global_[format_keys[i]].stream_id_local].load_order_index = i;
	}

	this->info_container_global_.clear();
	this->info_container_global_.resize(this->n_info_fields);
	this->info_container_local_.clear();
	this->info_container_local_.resize(block_footer.n_info_streams);
	this->info_container_loaded_.clear();

	for(U32 i = 0; i < block_footer.n_info_streams; ++i){
		// Set global -> local mapping -> loaded mapping
		this->info_container_global_[block_footer.info_offsets[i].data_header.global_key](-1, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
		// Set local -> global mapping -> loaded mapping
		this->info_container_local_[i](-1, i, block_footer.info_offsets[i].data_header.global_key, &block_footer.info_offsets[i]);
	}

	for(U32 i = 0; i < info_keys.size(); ++i){
		this->info_container_global_[info_keys[i]].load_order_index = i;
		this->info_container_local_[this->info_container_global_[info_keys[i]].stream_id_local].load_order_index = i;
	}

	return true;
}

}
}
