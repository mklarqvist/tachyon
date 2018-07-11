#ifndef CONTAINERS_VARIANT_READER_FILTERS_H_
#define CONTAINERS_VARIANT_READER_FILTERS_H_

#include "variant_reader_filters_tuple.h"
#include "variant_reader_objects.h"

namespace tachyon{

enum TACHYON_FILTER_FUNCTION{
	YON_FILTER_NUMBER_ALT_ALLELES,
	YON_FILTER_MIXED_PHASING,
	YON_FILTER_MIXED_PLOIDY,
	YON_FILTER_MISSING_GT,
	YON_FILTER_ALLELE_FREQUENCY,
	YON_FILTER_ALLELE_COUNT,
	YON_FILTER_UNIFORM_PHASE,
	YON_FILTER_KNOWN_NOVEL,
	YON_FILTER_REFERENCE_ALLELE,
	YON_FILTER_ALT_ALLELE,
	YON_FILTER_NAME
};

struct VariantReaderFilters{
public:
	typedef VariantReaderFilters self_type;
	typedef VariantReaderObjects objects_type;
	typedef VariantReaderFiltersTupleInterface value_type;
	typedef value_type&          reference;
	typedef const value_type&    const_reference;
	typedef value_type*          pointer;
	typedef const value_type*    const_pointer;
	typedef std::ptrdiff_t       difference_type;
	typedef std::size_t          size_type;
	typedef bool (self_type::*filter_function)(const_pointer pair, const objects_type& objects, const U32& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;

public:
	VariantReaderFilters() :
		n_filters_(0),
		n_capacity_(256),
		filter_data_(new pointer[this->n_capacity_]),
		require_genotypes(false),
		target_intervals(false)
	{

	}

	~VariantReaderFilters(){
		if(this->filter_data_ != nullptr){
			for(U32 i = 0; i < this->n_filters_; ++i)
				delete this->filter_data_[i];

			delete [] this->filter_data_;
		}
	}

	VariantReaderFilters(const VariantReaderFilters& other) = delete;

	template <class T>
	void add(TACHYON_FILTER_FUNCTION filter_function, const T& r_value, const TACHYON_COMPARATOR_TYPE& comparator){
		// Todo: currently if full then return: fix to resize and update
		if(this->size() + 1 == this->capacity())
			return;

		// Construct new filter function
		this->filter_data_[this->n_filters_++] = new VariantReaderFiltersTuple<T>(r_value, comparator);

		switch(filter_function){
		case(YON_FILTER_NUMBER_ALT_ALLELES):
			this->filters.push_back(&self_type::filterAlternativeAlleles);
			break;
		case(YON_FILTER_MIXED_PHASING):
			this->filters.push_back(&self_type::filterMixedPhasing);
			break;
		case(YON_FILTER_MIXED_PLOIDY):
			this->filters.push_back(&self_type::filterMixedPloidy);
			break;
		case(YON_FILTER_MISSING_GT):
			this->filters.push_back(&self_type::filterHasMissingGenotypes);
			break;
		case(YON_FILTER_ALLELE_FREQUENCY):
			this->filters.push_back(&self_type::filterAlleleFrequency);
			break;
		case(YON_FILTER_ALLELE_COUNT):
			this->filters.push_back(&self_type::filterAlleleCount);
			break;
		case(YON_FILTER_UNIFORM_PHASE):
			this->filters.push_back(&self_type::filterUniformMatchPhase);
			break;
		case(YON_FILTER_KNOWN_NOVEL):
			this->filters.push_back(&self_type::filterKnownNovel);
			break;
		case(YON_FILTER_REFERENCE_ALLELE):
			this->filters.push_back(&self_type::filterReferenceAllele);
			break;
		case(YON_FILTER_ALT_ALLELE):
			this->filters.push_back(&self_type::filterAlternativeAllele);
			break;
		case(YON_FILTER_NAME):
			this->filters.push_back(&self_type::filterName);
			break;
		}
	}

	inline const size_type& size(void) const{ return(this->n_filters_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Has mixed phasing
	inline bool filterMixedPhasing(const_pointer pair, const objects_type& objects, const U32& position) const{
		return(pair->applyFilter(objects.meta_container->at(position).isGTMixedPhasing()));
	}

	inline bool filterMixedPloidy(const_pointer pair, const objects_type& objects, const U32& position) const{
		return(pair->applyFilter((objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1])));
	}

	inline bool filterKnownNovel(const_pointer pair, const objects_type& objects, const U32& position) const{
		return(pair->applyFilter((U32)objects.meta_container->at(position).name.size()));
	}

	// GT data matches this
	inline bool filterUniformMatchPhase(const_pointer pair, const objects_type& objects, const U32& position) const
	{
		if(objects.meta_container->at(position).isGTMixedPhasing() == true) return false;
		return(pair->applyFilter(objects.meta_container->at(position).controller.gt_phase));
	}

	bool filterPloidy(const_pointer pair, const objects_type& objects, const U32& position) const;
	bool filterSampleList(const_pointer pair, const objects_type& objects, const U32& position) const;

	// BCFtools calculate this as the SUM of all ALT counts
	// We filter based on ANY ALT frequency OPERATOR the target frequency
	bool filterAlleleFrequency(const_pointer pair, const objects_type& objects, const U32& position) const{
		const std::vector<double> af = objects.genotype_summary->calculateAlleleFrequency(objects.meta_container->at(position));
		for(U32 i = 1; i < af.size(); ++i){
			if(pair->applyFilter(af[i]))
				return true;
		}
		return(false);
	}

	bool filterVariantClassification(const_pointer pair, const objects_type& object, const U32& position) const;
	bool filterUnseenAlternativeAlleles(const_pointer pair, const objects_type& object, const U32& position) const;
	bool filterFILTER(const_pointer pair, const objects_type& object, const U32& position) const;  // Filter by desired FILTER values
	bool filterINFO(const_pointer pair, const objects_type& object, const U32& position) const;    // custom filter. e.g. AC<1024

	inline bool filterAlternativeAlleles(const_pointer pair, const objects_type& object, const U32& position) const{
		// Remove one to total count as REF is counted here
		// Recast as signed integer to avoid possible underflowing issues
		return(pair->applyFilter(object.meta_container->at(position).getNumberAlleles() - 1));
	}

	inline bool filterAlleleCount(const_pointer pair, const objects_type& object, const U32& position) const{
		for(U32 i = 1; i < object.meta_container->at(position).n_alleles; ++i){
			if(pair->applyFilter(object.genotype_summary->vectorA_[2+i] + object.genotype_summary->vectorB_[2+i])){
				return true;
			}
		}
		return false;
	}

	inline bool filterHasMissingGenotypes(const_pointer pair, const objects_type& object, const U32& position) const{
		return(pair->applyFilter(object.genotype_summary->vectorA_[1]));
	}

	inline bool filterReferenceAllele(const_pointer pair, const objects_type& object, const U32& position) const{
		//std::cerr << object.meta->at(position).alleles[0].toString() << std::endl;
		return(pair->applyFilter(object.meta_container->at(position).alleles[0].toString()));
	}

	inline bool filterAlternativeAllele(const_pointer pair, const objects_type& object, const U32& position) const{
		for(U32 i = 1; i < object.meta_container->at(position).n_alleles; ++i){
			if(pair->applyFilter(object.meta_container->at(position).alleles[i].toString()))
				return true;
		}
		return false;
	}

	inline bool filterName(const_pointer pair, const objects_type& object, const U32& position) const{
		return(pair->applyFilter(object.meta_container->at(position).name));
	}

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool filter(const objects_type& objects, const U32 position) const{
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

	/**<
	 * Checks if any filter function require genotype data to be loaded and prepared
	 * @return Returns TRUE if genotype data is required or FALSE otherwise
	 */
	inline const bool doRequireGenotypes(void) const{ return(this->require_genotypes); }

public:
	size_type n_filters_;   // number of filters
	size_type n_capacity_;  // capacity
	bool require_genotypes; // Filtering require genotypes
	bool target_intervals;  // Filtering require intervals
	std::vector<filter_function> filters;
	value_type** filter_data_; // actual tuples stored here -> have to be double-pointer because of different payloads
};


}



#endif /* CONTAINERS_VARIANT_READER_FILTERS_H_ */
