#ifndef CORE_VARIANT_SITE_ANNOTATION_H_
#define CORE_VARIANT_SITE_ANNOTATION_H_

#include "genotype_summary.h"

namespace tachyon{
namespace containers{

struct VariantSiteAnnotation{
public:
	typedef VariantSiteAnnotation  self_type;
	typedef GenotypeSummary        genotype_summary_type;

public:
	VariantSiteAnnotation() :
		n_observed_alleles(0),
		n_eov(0),
		n_missing(0),
		genotype_summary(nullptr)
	{}

	~VariantSiteAnnotation(){}

	// Copy ctor
	VariantSiteAnnotation(const self_type& other) :
		n_observed_alleles(other.n_observed_alleles),
		n_eov(other.n_eov),
		n_missing(other.n_missing),
		hwe_p_values(other.hwe_p_values),
		allele_frequencies(other.allele_frequencies),
		strand_bias_fisher_phred_p(other.strand_bias_fisher_phred_p),
		genotype_summary(other.genotype_summary)
	{}

	// Todo: move ctor and assign move operator

public:
	U64 n_observed_alleles;
	U64 n_eov;
	U64 n_missing;
	std::vector<double>    hwe_p_values;
	std::vector<double>    allele_frequencies;
	std::vector<double>    strand_bias_fisher_phred_p;
	genotype_summary_type* genotype_summary;
};

}
}


#endif /* CORE_VARIANT_SITE_ANNOTATION_H_ */
