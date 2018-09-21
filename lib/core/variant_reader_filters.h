#ifndef CONTAINERS_VARIANT_READER_FILTERS_H_
#define CONTAINERS_VARIANT_READER_FILTERS_H_

#include "support/magic_constants.h"
#include "variant_record.h"
#include "variant_reader_filters_tuple.h"

namespace tachyon {

struct VariantReaderFilters {
public:
	typedef VariantReaderFilters self_type;
	typedef VariantReaderFiltersTupleInterface value_type;
	typedef value_type&          reference;
	typedef const value_type&    const_reference;
	typedef value_type*          pointer;
	typedef const value_type*    const_pointer;
	typedef std::ptrdiff_t       difference_type;
	typedef std::size_t          size_type;

	typedef bool (self_type::*filter_function)(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;

public:
	VariantReaderFilters();
	VariantReaderFilters(const VariantReaderFilters& other) = delete;
	VariantReaderFilters& operator=(const self_type& other) = delete;
	~VariantReaderFilters();

	template <class T>
	void Add(TACHYON_FILTER_FUNCTION filter_function,
	         const T& r_value,
	         const TACHYON_COMPARATOR_TYPE& comparator);

	// Capacity
	inline const size_type& size(void) const{ return(this->n_filters_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Has mixed phasing
	inline bool FilterMixedPhasing(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter(objects.controller.gt_has_mixed_phasing));
	}

	inline bool FilterKnownNovel(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter((uint32_t)objects.name.size()));
	}

	inline bool FilterQuality(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter(objects.qual));
	}

	// GT data matches this
	inline bool FilterUniformMatchPhase(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		if(objects.controller.gt_has_mixed_phasing == true) return false;
		return(pair->applyFilter(objects.controller.gt_has_mixed_phasing));
	}

	// BCFtools calculate this as the SUM of all ALT counts
	// We filter based on ANY ALT frequency OPERATOR the target frequency
	bool FilterAlleleFrequency(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool FilterUnseenAlternativeAlleles(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;

	inline bool FilterAlternativeAlleles(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		// Remove one to total count as REF is counted here
		// Recast as signed integer to avoid possible underflowing issues
		return(pair->applyFilter(objects.n_alleles - 1));
	}

	inline bool FilterAlleleCount(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i){
			if(pair->applyFilter((uint32_t)objects.gt_sum->d->ac[i])){
				return true;
			}
		}
		return false;
	}

	// Unused parameter position available in definition to allow a unified pointer definition
	inline bool FilterHasMissingGenotypes(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter((uint32_t)objects.gt_sum->d->ac[0]));
	}

	// Unused parameter position available in definition to allow a unified pointer definition
	inline bool FilterMixedPloidy(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter((uint32_t)objects.gt_sum->d->ac[1]));
	}

	inline bool FilterReferenceAllele(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		//std::cerr << objects.meta->at(position).alleles[0].toString() << std::endl;
		return(pair->applyFilter(objects.alleles[0].ToString()));
	}

	inline bool FilterAlternativeAllele(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		for(uint32_t i = 1; i < objects.n_alleles; ++i){
			if(pair->applyFilter(objects.alleles[i].ToString()))
				return true;
		}
		return false;
	}

	inline bool FilterName(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter(objects.name));
	}

	// Not implemented
	bool filterPloidy(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool filterSampleList(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool filterVariantClassification(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool filterFILTER(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;  // Filter by desired FILTER values
	bool filterINFO(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;    // custom filter. e.g. AC<1024

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool Filter(yon1_vnt_t& objects, const uint32_t position) const;

	/**<
	 * Checks if any filter function require genotype data to be loaded and prepared
	 * @return Returns TRUE if genotype data is required or FALSE otherwise
	 */
	inline bool HasRequireGenotypes(void) const{ return(this->require_genotypes); }

public:
	size_type n_filters_;   // number of filters
	size_type n_capacity_;  // capacity
	bool require_genotypes; // Filtering require genotypes
	bool target_intervals;  // Filtering require intervals
	std::vector<filter_function> filters;
	value_type** filter_data_; // actual tuples stored here -> have to be double-pointer because of different payloads
};


template <class T>
void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function,
                               const T& r_value,
                               const TACHYON_COMPARATOR_TYPE& comparator)
{
	// Todo: currently if full then return: fix to resize and update
	if(this->size() + 1 == this->capacity()){
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->filter_data_[this->n_filters_++] = new VariantReaderFiltersTuple<T>(r_value, comparator);

	switch(filter_function){
	case(YON_FILTER_NUMBER_ALT_ALLELES):
		this->filters.push_back(&self_type::FilterAlternativeAlleles);
		break;
	case(YON_FILTER_MIXED_PHASING):
		this->filters.push_back(&self_type::FilterMixedPhasing);
		break;
	case(YON_FILTER_MIXED_PLOIDY):
		this->filters.push_back(&self_type::FilterMixedPloidy);
		break;
	case(YON_FILTER_MISSING_GT):
		this->filters.push_back(&self_type::FilterHasMissingGenotypes);
		break;
	case(YON_FILTER_ALLELE_FREQUENCY):
		this->filters.push_back(&self_type::FilterAlleleFrequency);
		break;
	case(YON_FILTER_ALLELE_COUNT):
		this->filters.push_back(&self_type::FilterAlleleCount);
		break;
	case(YON_FILTER_UNIFORM_PHASE):
		this->filters.push_back(&self_type::FilterUniformMatchPhase);
		break;
	case(YON_FILTER_KNOWN_NOVEL):
		this->filters.push_back(&self_type::FilterKnownNovel);
		break;
	case(YON_FILTER_REFERENCE_ALLELE):
		this->filters.push_back(&self_type::FilterReferenceAllele);
		break;
	case(YON_FILTER_ALT_ALLELE):
		this->filters.push_back(&self_type::FilterAlternativeAllele);
		break;
	case(YON_FILTER_NAME):
		this->filters.push_back(&self_type::FilterName);
		break;
	case(YON_FILTER_UNSEEN_ALT):
		this->filters.push_back(&self_type::FilterUnseenAlternativeAlleles);
		break;
	case(YON_FILTER_QUALITY):
		this->filters.push_back(&self_type::FilterQuality);
		break;
	}
}


}



#endif /* CONTAINERS_VARIANT_READER_FILTERS_H_ */
