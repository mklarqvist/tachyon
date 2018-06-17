#ifndef CONTAINERS_VARIANT_READER_FILTERS_H_
#define CONTAINERS_VARIANT_READER_FILTERS_H_

#include "variant_reader_objects.h"

namespace tachyon{

template <class ValueClass>
struct VariantReaderFiltersTuple{
public:
	typedef VariantReaderFiltersTuple<ValueClass> self_type;
	typedef bool (self_type::*filter_function)(const ValueClass& target, const ValueClass& limit) const;

public:
	VariantReaderFiltersTuple() :
		filter(false),
		r_value(0),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value) :
		filter(true),
		r_value(r_value),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator) :
		filter(true),
		r_value(r_value),
		comparator(nullptr)
	{
		switch(comparator){
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		}
	}

	void operator()(const ValueClass& r_value){
		this->filter = true;
		this->r_value = r_value;
	}

	void operator()(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator){
		this->filter  = true;
		this->r_value = r_value;

		switch(comparator){
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		}
	}

	~VariantReaderFiltersTuple() = default;

	inline bool applyFilter(const ValueClass& l_value) const{ return((this->*comparator)(l_value, r_value)); }

	// Comparator functions
	inline bool __filterLesser(const ValueClass& target, const ValueClass& limit) const{return(target < limit);}
	inline bool __filterLesserEqual(const ValueClass& target, const ValueClass& limit) const{return(target <= limit);}
	inline bool __filterGreater(const ValueClass& target, const ValueClass& limit) const{return(target > limit);}
	inline bool __filterGreaterEqual(const ValueClass& target, const ValueClass& limit) const{return(target >= limit);}
	inline bool __filterEqual(const ValueClass& target, const ValueClass& limit) const{return(target == limit);}
	inline bool __filterNotEqual(const ValueClass& target, const ValueClass& limit) const{return(target != limit);}

public:
	bool            filter;
	ValueClass      r_value;
	filter_function comparator;
};

struct VariantReaderFilters{
public:
	typedef VariantReaderFilters self_type;
	typedef VariantReaderObjects objects_type;
	typedef bool (self_type::*filter_function)(const objects_type& objects, const U32& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;

	// example:
	// TACHYON_COMPARATOR_TYPE::YON_CMP_EQUAL

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
		return(this->filter_mixed_phase.applyFilter(objects.meta->at(position).isGTMixedPhasing()));
	}

	inline bool filterMixedPloidy(const objects_type& objects, const U32& position) const{
		//assert(objects.meta != nullptr);
		return(this->filter_mixed_ploidy.applyFilter(objects.meta->at(position).isAnyGTMixedPloidy()));
	}

	inline bool filterKnownNovel(const objects_type& objects, const U32& position) const{
		//assert(objects.meta != nullptr);
		return(this->filter_known_novel.applyFilter(objects.meta->at(position).name.size() == 0));
	}

	// GT data matches this
	inline bool filterUniformMatchPhase(const objects_type& objects,
												 const U32& position) const
	{
		//assert(objects.meta != nullptr);
		return(objects.meta->at(position).isGTMixedPhasing() == false &&
			   objects.meta->at(position).controller.gt_phase == this->filter_uniform_phase.r_value);
	}

	bool filterUncalled(const objects_type& objects, const U32& position) const;
	bool filterPloidy(const objects_type& objects, const U32& position) const;

	// Requires additional data
	bool filterSampleList(const objects_type& objects, const U32& position) const;


	// BCFtools calculate this as the SUM of all ALT counts
	// We filter based on ANY ALT frequency OPERATOR the target frequency
	bool filterAlleleFrequency(const objects_type& objects, const U32& position) const{
		const std::vector<double> af = objects.genotype_summary->calculateAlleleFrequency(objects.meta->at(position));
		for(U32 i = 1; i < af.size(); ++i){
			if(this->filter_af.applyFilter(af[i]))
				return true;
		}
		return(false);
	}

	bool filterVariantClassification(const objects_type& objects, const U32& position) const;
	bool filterUnseenAlternativeAlleles(const objects_type& objects, const U32& position) const;
	bool filterRegions(const objects_type& objects) const; // Filter by target intervals
	bool filterFILTER(const objects_type& objects) const;  // Filter by desired FILTER values
	bool filterINFO(const objects_type& objects) const;    // custom filter. e.g. AC<1024

	inline bool filterAlternativeAlleles(const objects_type& object, const U32& position) const{
		// Remove one to total count as REF is counted here
		// Recast as signed integer to avoid possible underflowing issues
		return(this->filter_n_alts.applyFilter(object.meta->at(position).getNumberAlleles() - 1));
	}

	inline bool filterAlleleCount(const objects_type& object, const U32& position) const{
		for(U32 i = 1; i < object.meta->at(position).n_alleles; ++i){
			if(this->filter_ac.applyFilter(object.genotype_summary->vectorA_[2+i] + object.genotype_summary->vectorB_[2+i])){
				return true;
			}
		}
		return false;
	}

	inline bool filterHasMissingGenotypes(const objects_type& object, const U32& position) const{
		return(this->filter_missing.applyFilter(object.genotype_summary->matrix_[0][0] == 0));
	}

	/**<
	 * Constructs the filter pointer vector given the fields that have been set
	 * @return Returns TRUE if passing construction or FALSE otherwise
	 */
	bool build(void){
		this->filters.clear();
		if(this->filter_n_alts.filter)        this->filters.push_back(&self_type::filterAlternativeAlleles);
		if(this->filter_mixed_phase.filter)   this->filters.push_back(&self_type::filterMixedPhasing);
		if(this->filter_mixed_ploidy.filter)  this->filters.push_back(&self_type::filterMixedPloidy);
		if(this->filter_missing.filter)       this->filters.push_back(&self_type::filterHasMissingGenotypes);
		if(this->filter_af.filter)            this->filters.push_back(&self_type::filterAlleleFrequency);
		if(this->filter_ac.filter)            this->filters.push_back(&self_type::filterAlleleCount);
		if(this->filter_uniform_phase.filter) this->filters.push_back(&self_type::filterUniformMatchPhase);
		if(this->filter_known_novel.filter)   this->filters.push_back(&self_type::filterKnownNovel);
		return true;
	}

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool filter(const objects_type& objects, const U32 position) const{
		if(this->require_genotypes)
			objects.genotypes->at(position).getSummary(*objects.genotype_summary);

		for(U32 i = 0 ; i < this->filters.size(); ++i){
			if((this->*(this->filters[i]))(objects, position) == false){
				return false;
			}
		}
		return true;
	}

public:
	bool require_genotypes;
	bool target_intervals;
	// std::vector<intervals> intervals;
	std::vector<filter_function>        filters;
	std::vector<family_filter_function> family_filters;

	VariantReaderFiltersTuple<bool>  filter_uniform_phase;
	VariantReaderFiltersTuple<SBYTE> filter_n_alts;
	VariantReaderFiltersTuple<bool>  filter_missing;
	VariantReaderFiltersTuple<bool>  filter_mixed_phase;
	VariantReaderFiltersTuple<bool>  filter_mixed_ploidy;
	VariantReaderFiltersTuple<bool>  filter_known_novel;
	VariantReaderFiltersTuple<float> filter_af;
	VariantReaderFiltersTuple<S32>   filter_ac;
};


}



#endif /* CONTAINERS_VARIANT_READER_FILTERS_H_ */
