#ifndef CONTAINERS_GENOTYPE_DATA_CONTAINER_INTERFACE_H_
#define CONTAINERS_GENOTYPE_DATA_CONTAINER_INTERFACE_H_

#include "../core/genotype_summary.h"
#include "../core/genotype_object.h"
#include "../core/ti_tv_object.h"
#include "../math/square_matrix.h"
#include "datacontainer.h"

namespace tachyon{
namespace containers{

#define YON_GT_RLE_ALLELE_A(PRIMITITVE, SHIFT, ADD) (((PRIMITITVE) & ((1 << (SHIFT)) - 1) << (ADD)) >> (ADD));
#define YON_GT_RLE_ALLELE_B(PRIMITIVE, SHIFT, ADD)  (((PRIMITIVE) & ((1 << (SHIFT)) - 1) << ((ADD)+(SHIFT))) >> ((ADD)+(SHIFT)));
#define YON_GT_RLE_LENGTH(PRIMITIVE, SHIFT, ADD) ((PRIMITIVE) >> (2*(SHIFT) + (ADD)))

class GenotypeContainerInterface{
private:
    typedef GenotypeContainerInterface       self_type;
    typedef std::size_t                      size_type;

protected:
    typedef tachyon::core::MetaEntry         meta_type;
    typedef tachyon::core::MetaHotController hot_controller_type;
    typedef tachyon::core::GTObject          gt_object;
    typedef GenotypeSum                      gt_summary;
    typedef math::SquareMatrix<double>       square_matrix_type;
    typedef algorithm::PermutationManager    permutation_type;
    typedef tachyon::core::TiTvObject        ti_tv_object_type;

public:
    GenotypeContainerInterface(void) :
    	n_entries(0),
		__data(nullptr),
		__meta(nullptr)
	{}

    GenotypeContainerInterface(const char* const data, const size_type& n_entries, const U32& n_bytes, const meta_type& meta) :
    	n_entries(n_entries),
		__data(new char[n_bytes]),
		__meta(&meta)
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

    /**<
     *
     * @param square_matrix
     * @return
     */
    virtual square_matrix_type& comparePairwise(square_matrix_type& square_matrix) const =0;

    /**<
     *
     * @return
     */
    virtual U32 getSum(void) const =0;

    // Summary statistics
    /**<
     *
     * @param gt_summary_object
     * @return
     */
    virtual gt_summary& updateSummary(gt_summary& gt_summary_object) const =0;

    /**<
     *
     * @return
     */
    virtual gt_summary getSummary(void) const =0;

    /**<
     *
     * @param gt_summary_object
     * @return
     */
    virtual gt_summary& getSummary(gt_summary& gt_summary_object) const =0;


    virtual std::vector<gt_object> getLiteralObjects(void) const =0;
    virtual std::vector<gt_object> getObjects(const U64& n_samples) const =0;
    virtual std::vector<gt_object> getObjects(const U64& n_samples, const permutation_type& ppa_manager) const =0;

    virtual void updateTransitionTransversions(std::vector<ti_tv_object_type>& objects) const =0;

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const meta_type& getMeta(void) const{ return(*this->__meta); }

protected:
    template <class Y>
    inline float comparatorSamplesDiploid(const Y& alleleA, const Y& ref_alleleA, const Y& alleleB, const Y& ref_alleleB) const{
    	// If allele A and allele B is identical in both samples but the alleles are different
		// e.g. 0|0 == 0|0
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
    const meta_type* __meta;
};

}
}

#endif /* CONTAINERS_GENOTYPE_DATA_CONTAINER_INTERFACE_H_ */
