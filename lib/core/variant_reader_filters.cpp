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
		for(U32 i = 0; i < this->n_filters_; ++i)
			delete this->filter_data_[i];

		delete [] this->filter_data_;
	}
}

bool VariantReaderFilters::filterAlleleFrequency(const_pointer pair, const objects_type& objects, const U32& position) const{
	const std::vector<double> af = objects.genotype_summary->calculateAlleleFrequency(objects.meta_container->at(position));
	for(U32 i = 1; i < af.size(); ++i){
		if(pair->applyFilter(af[i]))
			return true;
	}
	return(false);
}


bool VariantReaderFilters::filterUnseenAlternativeAlleles(const_pointer pair, const objects_type& object, const U32& position) const{
	for(U32 i = 0; i < object.meta_container->at(position).n_alleles; ++i){
		if(pair->applyFilter(object.genotype_summary->vectorA_[2+i] + object.genotype_summary->vectorB_[2+i] == 0))
			return true;
	}
	return false;
}

bool VariantReaderFilters::filter(const objects_type& objects, const U32 position) const{
	// Todo: construct genotype summary globally for this variant
	if(this->require_genotypes)
		objects.genotype_container->at(position).getSummary(*objects.genotype_summary);

	for(U32 i = 0 ; i < this->filters.size(); ++i){
		if((this->*(this->filters[i]))(this->filter_data_[i], objects, position) == false){
			return false;
		}
	}
	return true;
}


}
