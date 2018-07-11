#ifndef CONTAINERS_VARIANT_READER_FILTERS_H_
#define CONTAINERS_VARIANT_READER_FILTERS_H_

#include "variant_reader_filters_tuple.h"
#include "variant_reader_objects.h"

namespace tachyon{

struct VariantReaderFilters{
public:
	typedef VariantReaderFilters self_type;
	typedef VariantReaderObjects objects_type;

	typedef bool (self_type::*filter_function)(const objects_type& objects, const U32& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;



public:
	VariantReaderFilters() :
		require_genotypes(false),
		target_intervals(false)
	{

	}

	~VariantReaderFilters() = default;

	// Has mixed phasing
	inline bool filterMixedPhasing(const objects_type& objects, const U32& position) const{
		//assert(objects.meta != nullptr);
		return(this->filter_mixed_phase.applyFilter(objects.meta_container->at(position).isGTMixedPhasing()));
	}

	inline bool filterMixedPloidy(const objects_type& objects, const U32& position) const{
		//assert(objects.meta != nullptr);
		return(this->filter_mixed_ploidy.applyFilter((objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1])));
	}

	inline bool filterKnownNovel(const objects_type& objects, const U32& position) const{
		//assert(objects.meta != nullptr);
		return(this->filter_known_novel.applyFilter(objects.meta_container->at(position).name.size()));
	}

	// GT data matches this
	inline bool filterUniformMatchPhase(const objects_type& objects, const U32& position) const
	{
		//assert(objects.meta != nullptr);
		return(objects.meta_container->at(position).isGTMixedPhasing() == false &&
			   objects.meta_container->at(position).controller.gt_phase == this->filter_uniform_phase.r_value);
	}

	bool filterPloidy(const objects_type& objects, const U32& position) const;


	bool filterSampleList(const objects_type& objects, const U32& position) const;


	// BCFtools calculate this as the SUM of all ALT counts
	// We filter based on ANY ALT frequency OPERATOR the target frequency
	bool filterAlleleFrequency(const objects_type& objects, const U32& position) const{
		const std::vector<double> af = objects.genotype_summary->calculateAlleleFrequency(objects.meta_container->at(position));
		for(U32 i = 1; i < af.size(); ++i){
			if(this->filter_af.applyFilter(af[i]))
				return true;
		}
		return(false);
	}

	bool filterVariantClassification(const objects_type& object, const U32& position) const;
	bool filterUnseenAlternativeAlleles(const objects_type& object, const U32& position) const;
	bool filterRegions(const objects_type& object, const U32& position) const; // Filter by target intervals
	bool filterFILTER(const objects_type& object, const U32& position) const;  // Filter by desired FILTER values
	bool filterINFO(const objects_type& object, const U32& position) const;    // custom filter. e.g. AC<1024

	inline bool filterAlternativeAlleles(const objects_type& object, const U32& position) const{
		// Remove one to total count as REF is counted here
		// Recast as signed integer to avoid possible underflowing issues
		return(this->filter_n_alts.applyFilter(object.meta_container->at(position).getNumberAlleles() - 1));
	}

	inline bool filterAlleleCount(const objects_type& object, const U32& position) const{
		for(U32 i = 1; i < object.meta_container->at(position).n_alleles; ++i){
			if(this->filter_ac.applyFilter(object.genotype_summary->vectorA_[2+i] + object.genotype_summary->vectorB_[2+i])){
				return true;
			}
		}
		return false;
	}

	inline bool filterHasMissingGenotypes(const objects_type& object, const U32& position) const{
		return(this->filter_missing.applyFilter(object.genotype_summary->vectorA_[1]));
	}

	inline bool filterReferenceAllele(const objects_type& object, const U32& position) const{
		//std::cerr << object.meta->at(position).alleles[0].toString() << std::endl;
		return(this->filter_ref_allele.applyFilter(object.meta_container->at(position).alleles[0].toString()));
	}

	inline bool filterAlternativeAllele(const objects_type& object, const U32& position) const{
		for(U32 i = 1; i < object.meta_container->at(position).n_alleles; ++i){
			if(this->filter_alt_allele.applyFilter(object.meta_container->at(position).alleles[i].toString()))
				return true;
		}
		return false;
	}

	inline bool filterName(const objects_type& object, const U32& position) const{
		return(this->filter_name.applyFilter(object.meta_container->at(position).name));
	}

	/**<
	 * Constructs the filter pointer vector given the fields that have been set
	 */
	void build(void){
		this->filters.clear();
		if(this->filter_n_alts.filter)        this->filters.push_back(&self_type::filterAlternativeAlleles);
		if(this->filter_mixed_phase.filter)   this->filters.push_back(&self_type::filterMixedPhasing);
		if(this->filter_mixed_ploidy.filter)  this->filters.push_back(&self_type::filterMixedPloidy);
		if(this->filter_missing.filter)       this->filters.push_back(&self_type::filterHasMissingGenotypes);
		if(this->filter_af.filter)            this->filters.push_back(&self_type::filterAlleleFrequency);
		if(this->filter_ac.filter)            this->filters.push_back(&self_type::filterAlleleCount);
		if(this->filter_uniform_phase.filter) this->filters.push_back(&self_type::filterUniformMatchPhase);
		if(this->filter_known_novel.filter)   this->filters.push_back(&self_type::filterKnownNovel);
		if(this->filter_ref_allele.filter)    this->filters.push_back(&self_type::filterReferenceAllele);
		if(this->filter_alt_allele.filter)    this->filters.push_back(&self_type::filterAlternativeAllele);
		if(this->filter_name.filter)          this->filters.push_back(&self_type::filterName);
	}

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool filter(const objects_type& objects, const U32 position) const{
		if(this->require_genotypes)
			objects.genotype_container->at(position).getSummary(*objects.genotype_summary);

		for(U32 i = 0 ; i < this->filters.size(); ++i){
			if((this->*(this->filters[i]))(objects, position) == false){
				return false;
			}
		}
		return true;
	}

	/**<
	 * Checks if any filter function require genotype data to be loaded and prepared
	 * @return Returns TRUE if genotype data is required or FALSE otherwise
	 */
	inline const bool doRequireGenotypes(void) const{ return(this->require_genotypes); }

public:
	bool require_genotypes; // Filtering require genotypes
	bool target_intervals;  // Filtering require intervals
	// std::vector<intervals> intervals;
	std::vector<filter_function>        filters;
	std::vector<family_filter_function> family_filters;

	VariantReaderFiltersTuple<bool>  filter_uniform_phase;
	VariantReaderFiltersTuple<SBYTE> filter_n_alts;
	VariantReaderFiltersTuple<S64>   filter_missing;
	VariantReaderFiltersTuple<bool>  filter_mixed_phase;
	VariantReaderFiltersTuple<bool>  filter_mixed_ploidy;
	VariantReaderFiltersTuple<bool>  filter_known_novel;
	VariantReaderFiltersTuple<float> filter_af;
	VariantReaderFiltersTuple<S32>   filter_ac;

	VariantReaderFiltersTuple<std::string>   filter_ref_allele;
	VariantReaderFiltersTuple<std::string>   filter_alt_allele;
	VariantReaderFiltersTuple<std::string>   filter_name;
};


}



#endif /* CONTAINERS_VARIANT_READER_FILTERS_H_ */
