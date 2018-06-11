#ifndef CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_
#define CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_

#include "../algorithm/permutation/permutation_manager.h"
#include "../core/genotype_object.h"
#include "../core/genotype_summary.h"
#include "../core/ts_tv_object.h"
#include "../math/square_matrix.h"
#include "datacontainer.h"

namespace tachyon{
namespace containers{

// Todo: value type should be GT object
class GenotypeContainerInterface{
private:
    typedef GenotypeContainerInterface    self_type;
    typedef std::size_t                   size_type;

protected:
    typedef core::MetaEntry               meta_type;
    typedef core::VariantController       hot_controller_type;
    typedef core::GTObject                gt_object;
    typedef GenotypeSummary               gt_summary;
    typedef math::SquareMatrix<double>    square_matrix_type;
    typedef algorithm::PermutationManager permutation_type;
    typedef core::TsTvObject              ts_tv_object_type;

    // Function pointers
	typedef float (self_type::*matrix_comparator)(const BYTE& alleleA, const BYTE& ref_alleleA, const BYTE& alleleB, const BYTE& ref_alleleB);

public:
    GenotypeContainerInterface(void) :
    	n_entries(0),
		__data(nullptr)
	{}

    GenotypeContainerInterface(const char* const data, const size_type& n_entries, const U32& n_bytes, const meta_type& meta) :
    	n_entries(n_entries),
		__data(new char[n_bytes]),
		__meta(meta)
    {
    	memcpy(this->__data, data, n_bytes);
    }

    virtual ~GenotypeContainerInterface(){ delete [] this->__data; }

    // GT-specific functionality
    //virtual void getGTSummary(void) const =0;
    //virtual void getGTSummaryGroups(void) const =0;
    //virtual void std::vector<float> getAlleleFrequency(void) =0;
    //virtual void std::vector<bool> getSamplesMissingness(void) =0;
    //virtual void std::vector<U32> getSamplesPloidy(void) =0;
    //virtual void std::vector<sample_summary> getSamplesSummary(void) =0;

    virtual square_matrix_type& comparePairwise(square_matrix_type& square_matrix) const =0;
    virtual U32 getSum(void) const =0;

    // Summary statistics
    virtual gt_summary& updateSummary(gt_summary& gt_summary_object) const =0;
    virtual gt_summary getSummary(void) const =0;
    virtual gt_summary& getSummary(gt_summary& gt_summary_object) const =0;


    virtual std::vector<gt_object> getLiteralObjects(void) const =0;
    virtual std::vector<gt_object> getObjects(const U64& n_samples) const =0;
    virtual std::vector<gt_object> getObjects(const U64& n_samples, const permutation_type& ppa_manager) const =0;
	virtual void getObjects(const U64& n_samples, std::vector<gt_object>& objects) const =0;
	virtual void getObjects(const U64& n_samples, std::vector<gt_object>& objects, const permutation_type& ppa_manager) const =0;
	virtual void getLiteralObjects(std::vector<gt_object>& objects) const =0;
    virtual void getTsTv(std::vector<ts_tv_object_type>& objects) const =0;

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const meta_type& getMeta(void) const{ return(this->__meta); }

protected:
    inline float comparatorSamplesDiploid(const BYTE& alleleA, const BYTE& ref_alleleA, const BYTE& alleleB, const BYTE& ref_alleleB) const{
    	// If allele A and allele B is identical in both samples but the alleles are different
		// e.g. 0|1 == 1|0
    	if((alleleA == ref_alleleA && alleleB == ref_alleleB) || (alleleA == ref_alleleB && alleleB == ref_alleleA)){ // identical but heterozygote
			return(1);
		}
    	// If allele A and allele B are identical in both samples
		// e.g. 0|0 == 0|0
    	else if(alleleA == ref_alleleA && alleleB == ref_alleleB){ // identical
			return(2);
		}
		// If neither allele A nor B are the same
		// e.g. 0|1 == 1|0
		// e.g. 1|1 == 0|0
		else if(!(alleleA == ref_alleleA || alleleB == ref_alleleB)){
			return(0);
		}
		// Otherwise one allele matches
		// e.g. 0|0 == 0|1
		// e.g. 1|0 == 1|1
		else {
			return(1);
		}
    }

protected:
    size_type        n_entries;
    char*            __data;
    const meta_type  __meta;
};

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_ */
