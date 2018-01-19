#include "MetaIterator.h"

namespace Tachyon{
namespace Iterator{

MetaIterator::MetaIterator() :
	current_position(0),
	n_entries(0),
	n_entries_unfiltered(0),
	gt_filter_option(Core::YON_GT_KEEP_ALL),
	container_hot(nullptr),
	container_cold(nullptr),
	info_id_container(nullptr),
	filter_id_container(nullptr),
	format_id_container(nullptr)
{

}

MetaIterator::MetaIterator(const container_type& container) :
	current_position(0),
	n_entries(0),
	n_entries_unfiltered(0),
	gt_filter_option(Core::YON_GT_KEEP_ALL),
	hot_iterator(container),
	container_hot(&container),
	container_cold(nullptr),
	info_id_container(nullptr),
	filter_id_container(nullptr),
	format_id_container(nullptr)

{
	this->setup(container);
}

MetaIterator::MetaIterator(const container_type& container, const Core::TACHYON_GT_TYPE gt_filter_type) :
	current_position(0),
	n_entries(0),
	n_entries_unfiltered(0),
	gt_filter_option(gt_filter_type),
	hot_iterator(container),
	container_hot(&container),
	container_cold(nullptr),
	info_id_container(nullptr),
	filter_id_container(nullptr),
	format_id_container(nullptr)

{
	this->setup(container);
}

MetaIterator::MetaIterator(const container_type& container_hot, const container_type& container_cold) :
	current_position(0),
	n_entries(0),
	n_entries_unfiltered(0),
	gt_filter_option(Core::YON_GT_KEEP_ALL),
	hot_iterator(container_hot),
	cold_iterator(container_cold, hot_iterator.size()),
	container_hot(&container_hot),
	container_cold(&container_cold),
	info_id_container(nullptr),
	filter_id_container(nullptr),
	format_id_container(nullptr)
{
	this->setup(container_hot, container_cold);
}

MetaIterator::MetaIterator(const container_type& container_hot, const container_type& container_cold, const Core::TACHYON_GT_TYPE gt_filter_type) :
	current_position(0),
	n_entries(0),
	n_entries_unfiltered(0),
	gt_filter_option(gt_filter_type),
	hot_iterator(container_hot),
	cold_iterator(container_cold, hot_iterator.size()),
	container_hot(&container_hot),
	container_cold(&container_cold),
	info_id_container(nullptr),
	filter_id_container(nullptr),
	format_id_container(nullptr)
{
	this->setup(container_hot, container_cold);
}

MetaIterator::~MetaIterator(){
	this->clearPrevious();
}

bool MetaIterator::setup(const container_type& hot_meta_container){
	if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
		return false;

	this->container_hot = &hot_meta_container;
	if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
		this->hot_iterator.setup(hot_meta_container);
	else
		this->hot_iterator.setup(hot_meta_container, this->gt_filter_option);

	this->clearPrevious();

	assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
	this->n_entries = this->hot_iterator.size();
	this->n_entries_unfiltered = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);

	this->entries.resize(this->n_entries);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->entries[i] = entry_type(this->hot_iterator[i]);

	return true;
}

bool MetaIterator::setup(const container_type& hot_meta_container, const container_type& cold_meta_container){
	if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
		return false;

	this->container_hot  = &hot_meta_container;
	this->container_cold = &cold_meta_container;
	if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
		this->hot_iterator.setup(hot_meta_container);
	else
		this->hot_iterator.setup(hot_meta_container, this->gt_filter_option);

	this->clearPrevious();

	assert((this->container_hot->buffer_data_uncompressed.pointer % sizeof(hot_type)) == 0);
	this->n_entries_unfiltered = this->container_hot->buffer_data_uncompressed.pointer / sizeof(hot_type);
	this->n_entries = this->hot_iterator.size();
	if(this->gt_filter_option == Core::YON_GT_KEEP_ALL)
		this->cold_iterator.setup(cold_meta_container, this->n_entries_unfiltered);
	else
		this->cold_iterator.setup(cold_meta_container, this->n_entries_unfiltered, this->gt_filter_option);

	this->entries.resize(this->n_entries);
	for(U32 i = 0; i < this->n_entries; ++i)
		this->entries[i] = entry_type(this->hot_iterator[i], this->cold_iterator[i]);

	return true;
}

bool MetaIterator::setInfoIDContainer(const container_type& info_id_container){
	if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
		return false;

	if(this->n_entries <= 0) return false;
	if(this->n_entries_unfiltered <=0) return false;

	this->info_id_container = &info_id_container;
	this->info_id_iterator.setup(info_id_container);
	const void* p;
	if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
		const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
		U32 j = 0;
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->info_id_iterator.getDataIterator()->currentPointer(p);
			if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
				this->entries[j++].getInfoPatternID() = *(S32*)p;
			++this->info_id_iterator;
		}
	} else {
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->info_id_iterator.getDataIterator()->currentPointer(p);
			this->entries[i].getInfoPatternID() = *(S32*)p;
			++this->info_id_iterator;
		}
	}
	return true;
}

bool MetaIterator::setFilterIDContainer(const container_type& filter_id_container){
	if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
		return false;

	if(this->n_entries == 0) return false;
	if(this->n_entries_unfiltered <=0) return false;

	this->filter_id_container = &filter_id_container;
	this->filter_id_iterator.setup(filter_id_container);
	const void* p;
	if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
		const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
		U32 j = 0;
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->filter_id_iterator.getDataIterator()->currentPointer(p);
			if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
				this->entries[j++].getFilterPatternID() = *(S32*)p;
			++this->filter_id_iterator;
		}
	} else {
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->filter_id_iterator.getDataIterator()->currentPointer(p);
			this->entries[i].getFilterPatternID() = *(S32*)p;
			++this->filter_id_iterator;
		}
	}
	return true;
}

bool MetaIterator::setFormatIDContainer(const container_type& format_id_container){
	if(this->gt_filter_option == Core::YON_GT_UNKNOWN)
		return false;

	if(this->n_entries == 0) return false;
	if(this->n_entries_unfiltered <=0) return false;

	this->format_id_container = &format_id_container;
	this->format_id_iterator.setup(format_id_container);
	const void* p;
	if(this->gt_filter_option != Core::YON_GT_KEEP_ALL){
		const hot_type* const hot_entries_all = reinterpret_cast<const hot_type* const>(this->container_hot->buffer_data_uncompressed.data);
		U32 j = 0;
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->format_id_iterator.getDataIterator()->currentPointer(p);
			if(hot_entries_all[i].getGenotypeType() == this->gt_filter_option)
				this->entries[j++].getFormatPatternID() = *(S32*)p;
			++this->format_id_iterator;
		}
	} else {
		for(U32 i = 0; i < this->n_entries_unfiltered; ++i){
			this->format_id_iterator.getDataIterator()->currentPointer(p);
			this->entries[i].getFormatPatternID() = *(S32*)p;
			++this->format_id_iterator;
		}
	}
	return true;
}

void MetaIterator::reset(void){
	this->current_position = 0;
	this->hot_iterator.reset();
	this->cold_iterator.reset();
	this->info_id_iterator.reset();
	this->filter_id_iterator.reset();
	this->format_id_iterator.reset();
}

}
}
