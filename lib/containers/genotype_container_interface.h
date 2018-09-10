#ifndef CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_
#define CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_

#include "data_container.h"
#include "core/genotypes.h"

namespace tachyon {
namespace containers {

class GenotypeContainerInterface {
public:
    typedef GenotypeContainerInterface    self_type;
    typedef std::size_t                   size_type;
    typedef core::MetaEntry               meta_type;
    typedef yon_gt_summary                gt_summary;

    // Function pointers
	typedef float (self_type::*matrix_comparator)(const uint8_t& alleleA, const uint8_t& ref_alleleA, const uint8_t& alleleB, const uint8_t& ref_alleleB);

public:
    GenotypeContainerInterface(void) :
    	n_entries(0),
		__data(nullptr),
		gt_missing(false),
		gt_mixed_phase(false),
		gt_global_phase(false),
		n_base_ploidy(0),
		n_alleles(0)
	{}

    GenotypeContainerInterface(const char* const data, const size_type n_entries, const uint32_t& n_bytes, const meta_type& meta) :
    	n_entries(n_entries),
		__data(new uint8_t[n_bytes]),
		gt_missing(meta.controller.gt_anyMissing),
		gt_mixed_phase(meta.controller.gt_mixed_phasing),
		gt_global_phase(meta.controller.gt_phase),
		n_base_ploidy(meta.n_base_ploidy),
		n_alleles(meta.n_alleles)
    {
    	memcpy(this->__data, data, n_bytes);
    }

    GenotypeContainerInterface(const char* const data,
    		const size_type& n_entries,
			const uint32_t& n_bytes,
			const uint16_t n_als,
			const uint8_t  n_ploidy,
			const uint16_t controller) :
		n_entries(n_entries),
		__data(new uint8_t[n_bytes])
	{
    	const yon_vnt_cnt* cont = reinterpret_cast<const yon_vnt_cnt*>(&controller);
		gt_missing = cont->gt_has_missing;
		gt_mixed_phase = cont->gt_has_mixed_phasing;
		gt_global_phase = cont->gt_phase_uniform;
		n_base_ploidy = n_ploidy;
		n_alleles = n_als;
		memcpy(this->__data, data, n_bytes);
	}

    virtual ~GenotypeContainerInterface(){ delete [] this->__data; }

    virtual yon_gt* GetObjects(const uint32_t n_samples) =0;
    virtual yon_gt* GetObjects(yon_gt_ppa& ppa) =0;

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

protected:
    size_type        n_entries;
    uint8_t*         __data;

    bool gt_missing;
    bool gt_mixed_phase;
    bool gt_global_phase;
    uint8_t n_base_ploidy;
    uint16_t n_alleles;
};

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_ */
