#include "variant_reader_filters.h"

namespace tachyon{

VariantReaderFilters::VariantReaderFilters() :
	n_filters_(0),
	n_capacity_(256),
	require_genotypes(false),
	target_intervals(false),
	filter_data_(new pointer[this->n_capacity_])
{

}

VariantReaderFilters::~VariantReaderFilters(){
	if(this->filter_data_ != nullptr){
		for(uint32_t i = 0; i < this->n_filters_; ++i)
			delete this->filter_data_[i];

		delete [] this->filter_data_;
	}
}

bool VariantReaderFilters::filterAlleleFrequency(const_pointer pair, const yon1_t& objects, const uint32_t& position) const{
	for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i){
		if(pair->applyFilter(objects.gt_sum->d->af[i]))
			return true;
	}
	return(false);
}


bool VariantReaderFilters::filterUnseenAlternativeAlleles(const_pointer pair, const yon1_t& objects, const uint32_t& position) const{
	for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i){
		if(pair->applyFilter(objects.gt_sum->d->ac[i] + objects.gt_sum->d->ac[i] == 0))
			return true;
	}
	return false;
}

bool VariantReaderFilters::filter(yon1_t& objects, const uint32_t position) const{
	if(this->require_genotypes)
		objects.EvaluateSummary(true);

	for(uint32_t i = 0 ; i < this->filters.size(); ++i){
		if((this->*(this->filters[i]))(this->filter_data_[i], objects, position) == false){
			return false;
		}
	}
	return true;
}


}
