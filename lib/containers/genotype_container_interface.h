#ifndef CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_
#define CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_

#include "data_container.h"
#include "core/genotypes.h"

namespace tachyon {
namespace containers {

class GenotypeContainerInterface{
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
		__data(nullptr)
	{}

    GenotypeContainerInterface(const char* const data, const size_type& n_entries, const uint32_t& n_bytes, const meta_type& meta) :
    	n_entries(n_entries),
		__data(new uint8_t[n_bytes]),
		__meta(meta)
    {
    	memcpy(this->__data, data, n_bytes);
    }

    virtual ~GenotypeContainerInterface(){ delete [] this->__data; }

    virtual yon_gt* GetObjects(const uint32_t n_samples) =0;
    virtual yon_gt* GetObjects(yon_gt_ppa& ppa) =0;

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const meta_type& GetMeta(void) const{ return(this->__meta); }

protected:
    size_type        n_entries;
    uint8_t*         __data;
    const meta_type  __meta;
};

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_INTERFACE_H_ */
