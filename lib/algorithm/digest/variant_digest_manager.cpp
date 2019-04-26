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
	for (uint32_t i = 0; i < this->n_entries_info_; ++i) this->__entries_info[i] = other.__entries_info[i];
	for (uint32_t i = 0; i < this->n_entries_format_; ++i) this->__entries_format[i] = other.__entries_format[i];
}

VariantDigestManager::~VariantDigestManager() {
	delete [] this->__entries_info;
	delete [] this->__entries_format;
}

void VariantDigestManager::Finalize(void) {
	parent_type::finalize();
	for (uint32_t i = 0; i < this->n_capacity_info_; ++i)   this->atINFO(i).finalize();
	for (uint32_t i = 0; i < this->n_capacity_format ; ++i) this->atFORMAT(i).finalize();
}

void VariantDigestManager::operator+=(const variant_block_type& block) {
	for (uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)
		this->at(i) += block.base_containers[i];

	for (uint32_t i = 0; i < block.footer.n_info_streams; ++i)   this->__entries_info[block.footer.info_offsets[i].data_header.global_key] += block.info_containers[i];
	for (uint32_t i = 0; i < block.footer.n_format_streams; ++i) this->__entries_format[block.footer.format_offsets[i].data_header.global_key] += block.format_containers[i];
}

}
}
