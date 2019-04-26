#include "variant_reader_filters_tuple.h"
#include "variant_reader_filters.h"

namespace tachyon {

class VariantReaderFilters::VariantReaderFiltersImpl {
public:
	typedef VariantReaderFiltersImpl self_type;
	typedef VariantReaderFiltersTupleInterface value_type;
	typedef value_type&       reference;
	typedef const value_type& const_reference;
	typedef value_type*       pointer;
	typedef const value_type* const_pointer;
	typedef std::ptrdiff_t    difference_type;
	typedef std::size_t       size_type;

	typedef bool (self_type::*filter_function)(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;

public:
	VariantReaderFiltersImpl() :
		n_filters_(0),
		n_capacity_(256),
		require_genotypes(false),
		target_intervals(false),
		filter_data_(new pointer[this->n_capacity_])
	{

	}

	~VariantReaderFiltersImpl() {
		if (this->filter_data_ != nullptr) {
			for(uint32_t i = 0; i < this->n_filters_; ++i)
				delete this->filter_data_[i];

			delete [] this->filter_data_;
		}
	}

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
		if (objects.controller.gt_has_mixed_phasing == true) return false;
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
		for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i) {
			if (pair->applyFilter((uint32_t)objects.gt_sum->d->ac[i])) {
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
		for(uint32_t i = 1; i < objects.n_alleles; ++i) {
			if (pair->applyFilter(objects.alleles[i].ToString()))
				return true;
		}
		return false;
	}

	inline bool FilterName(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
		return(pair->applyFilter(objects.name));
	}

	// Not implemented
	bool FilterPloidy(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool FilterSampleList(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool FilterVariantClassification(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;
	bool FilterFILTER(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;  // Filter by desired FILTER values
	bool FilterINFO(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const;    // custom filter. e.g. AC<1024

public:
	size_type n_filters_;   // number of filters
	size_type n_capacity_;  // capacity
	bool require_genotypes; // Filtering require genotypes
	bool target_intervals;  // Filtering require intervals

	std::vector<filter_function> filters;
	value_type** filter_data_; // actual tuples stored here -> have to be double-pointer because of different payloads

};

VariantReaderFilters::VariantReaderFilters() :
	mImpl(new VariantReaderFiltersImpl)
{

}

VariantReaderFilters::~VariantReaderFilters() {

}

// Capacity
const size_t& VariantReaderFilters::size(void) const { return(this->mImpl->n_filters_); }
const size_t& VariantReaderFilters::capacity(void) const { return(this->mImpl->n_capacity_); }

bool VariantReaderFilters::HasRequireGenotypes(void) const { return(this->mImpl->require_genotypes); }

void VariantReaderFilters::SetRequireGenotypes(bool set) { this->mImpl->require_genotypes = set; }

void VariantReaderFilters::AddWrapper(TACHYON_FILTER_FUNCTION filter_function) {
	switch(filter_function) {
	case(YON_FILTER_NUMBER_ALT_ALLELES):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterAlternativeAlleles);
		break;
	case(YON_FILTER_MIXED_PHASING):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterMixedPhasing);
		break;
	case(YON_FILTER_MIXED_PLOIDY):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterMixedPloidy);
		break;
	case(YON_FILTER_MISSING_GT):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterHasMissingGenotypes);
		break;
	case(YON_FILTER_ALLELE_FREQUENCY):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterAlleleFrequency);
		break;
	case(YON_FILTER_ALLELE_COUNT):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterAlleleCount);
		break;
	case(YON_FILTER_UNIFORM_PHASE):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterUniformMatchPhase);
		break;
	case(YON_FILTER_KNOWN_NOVEL):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterKnownNovel);
		break;
	case(YON_FILTER_REFERENCE_ALLELE):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterReferenceAllele);
		break;
	case(YON_FILTER_ALT_ALLELE):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterAlternativeAllele);
		break;
	case(YON_FILTER_NAME):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterName);
		break;
	case(YON_FILTER_UNSEEN_ALT):
		this->SetRequireGenotypes(true);
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterUnseenAlternativeAlleles);
		break;
	case(YON_FILTER_QUALITY):
		this->mImpl->filters.push_back(&VariantReaderFiltersImpl::FilterQuality);
		break;
	}
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const char& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<char>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const int8_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<int8_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const int16_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<int16_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const int32_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<int32_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const int64_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<int64_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const uint8_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<uint8_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const uint16_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<uint16_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const uint32_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<uint32_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const uint64_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<uint64_t>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const std::string& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<std::string>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const float& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<float>(r_value, comparator);
}

void VariantReaderFilters::Add(TACHYON_FILTER_FUNCTION filter_function, const double& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
	// Todo: currently if full then return: fix to resize and update
	if (this->size() + 1 == this->capacity()) {
		std::cerr << utility::timestamp("ERROR","KNOWN-BUG") << "Filter array is full..." << std::endl;
		return;
	}

	// Construct new filter function
	this->AddWrapper(filter_function);
	this->mImpl->filter_data_[this->mImpl->n_filters_++] = new VariantReaderFiltersTuple<double>(r_value, comparator);
}

bool VariantReaderFilters::VariantReaderFiltersImpl::FilterAlleleFrequency(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
	for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i) {
		if (pair->applyFilter(objects.gt_sum->d->af[i]))
			return true;
	}
	return(false);
}


bool VariantReaderFilters::VariantReaderFiltersImpl::FilterUnseenAlternativeAlleles(const_pointer pair, const yon1_vnt_t& objects, const uint32_t& position) const {
	for(uint32_t i = 3; i < objects.gt_sum->d->n_ac_af; ++i) {
		if (pair->applyFilter(objects.gt_sum->d->ac[i] + objects.gt_sum->d->ac[i] == 0))
			return true;
	}
	return false;
}

bool VariantReaderFilters::Filter(yon1_vnt_t& objects, const uint32_t position) const {
	if (this->mImpl->require_genotypes)
		objects.EvaluateSummary(true);

	// Iterate over the vector of filter function pointers and
	// invoke them in order.
	for(uint32_t i = 0 ; i < this->mImpl->filters.size(); ++i) {
		if ((this->mImpl.get()->*(this->mImpl->filters[i]))(this->mImpl->filter_data_[i], objects, position) == false) {
			return false;
		}
	}
	return true;
}


}
