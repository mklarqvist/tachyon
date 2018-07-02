#include "variant_digest_manager.h"

namespace tachyon{
namespace algorithm{

VariantDigestManager::VariantDigestManager() :
	parent_type(100),
	n_entries_info_(0),
	n_entries_format_(0),
	n_capacity_info_(100),
	n_capacity_format(100),
	__entries_info(new value_type[this->n_entries_info_]),
	__entries_format(new value_type[this->n_entries_format_])
{
}

VariantDigestManager::VariantDigestManager(const size_type base_capacity) :
	parent_type(base_capacity),
	n_entries_info_(base_capacity),
	n_entries_format_(base_capacity),
	n_capacity_info_(base_capacity),
	n_capacity_format(base_capacity),
	__entries_info(new value_type[this->n_entries_info_]),
	__entries_format(new value_type[this->n_entries_format_])
{
}

VariantDigestManager::VariantDigestManager(const size_type base_capacity, const size_type capacity_info, const size_type capacity_format) :
	parent_type(base_capacity),
	n_entries_info_(capacity_info),
	n_entries_format_(capacity_format),
	n_capacity_info_(capacity_info),
	n_capacity_format(capacity_format),
	__entries_info(new value_type[this->n_entries_info_]),
	__entries_format(new value_type[this->n_entries_format_])
{
}

VariantDigestManager::VariantDigestManager(const self_type& other) :
	parent_type(other),
	n_entries_info_(other.n_entries_info_),
	n_entries_format_(other.n_entries_format_),
	n_capacity_info_(other.n_capacity_info_),
	n_capacity_format(other.n_capacity_format),
	__entries_info(new value_type[this->n_entries_info_]),
	__entries_format(new value_type[this->n_entries_format_])
{
	for(U32 i = 0; i < this->n_entries_info_; ++i) this->__entries_info[i] = other.__entries_info[i];
	for(U32 i = 0; i < this->n_entries_format_; ++i) this->__entries_format[i] = other.__entries_format[i];
}

VariantDigestManager::~VariantDigestManager(){
	delete [] this->__entries_info;
	delete [] this->__entries_format;
}

void VariantDigestManager::finalize(void){
	parent_type::finalize();
	for(U32 i = 0; i < this->n_capacity_info_; ++i) this->atINFO(i).finalize();
	for(U32 i = 0; i < this->n_capacity_format ; ++i) this->atFORMAT(i).finalize();
}

void VariantDigestManager::operator+=(const variant_block_type& block){
	this->at(1)  += block.meta_contig_container;
	this->at(2)  += block.meta_positions_container;
	this->at(3)  += block.meta_names_container;
	this->at(4)  += block.meta_refalt_container;
	this->at(5)  += block.meta_controller_container;
	this->at(6)  += block.meta_quality_container;
	this->at(7)  += block.meta_names_container;
	this->at(8)  += block.meta_alleles_container;
	this->at(9)  += block.meta_info_map_ids;
	this->at(10) += block.meta_format_map_ids;
	this->at(11) += block.meta_filter_map_ids;
	this->at(12) += block.gt_support_data_container;
	this->at(13) += block.gt_rle8_container;
	this->at(14) += block.gt_rle16_container;
	this->at(15) += block.gt_rle32_container;
	this->at(16) += block.gt_rle64_container;
	this->at(17) += block.gt_simple8_container;
	this->at(18) += block.gt_simple16_container;
	this->at(19) += block.gt_simple32_container;
	this->at(20) += block.gt_simple64_container;

	for(U32 i = 0; i < block.footer.n_info_streams; ++i) this->__entries_info[block.footer.info_offsets[i].data_header.global_key] += block.info_containers[i];
	for(U32 i = 0; i < block.footer.n_format_streams; ++i) this->__entries_format[block.footer.format_offsets[i].data_header.global_key] += block.format_containers[i];
}

}
}
