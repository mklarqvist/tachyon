#ifndef CONTAINERS_GENOTYPE_CONTAINER_H_
#define CONTAINERS_GENOTYPE_CONTAINER_H_

#include "genotype_container_diploid_bcf.h"
#include "genotype_container_diploid_rle.h"
#include "genotype_container_diploid_simple.h"
#include "genotype_container_nploid.h"
#include "variant_block.h"

namespace tachyon{
namespace containers{

// forward declaration for getters
template <class T> class GenotypeContainerDiploidRLE;
template <class T> class GenotypeContainerDiploidSimple;


class GenotypeContainer {
public:
    typedef GenotypeContainer          self_type;
    typedef GenotypeContainerInterface value_type;
    typedef value_type&                reference;
    typedef const value_type&          const_reference;
    typedef value_type*                pointer;
    typedef const value_type*          const_pointer;
    typedef std::ptrdiff_t             difference_type;
    typedef std::size_t                size_type;
    typedef io::BasicBuffer            buffer_type;
    typedef yon_gt_summary             gt_summary_type;
    typedef VariantBlock               block_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    GenotypeContainer();
    ~GenotypeContainer();

    // Capacity
	inline bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Element access
	inline pointer         data(void){ return(this->__iterators); }
	inline const_pointer   data(void) const{ return(this->__iterators); }
	inline reference       operator[](const uint32_t& position){ return(this->__iterators[position]); }
	inline const_reference operator[](const uint32_t& position) const{ return(this->__iterators[position]); }
	inline reference       at(const uint32_t& position){ return(this->__iterators[position]); }
	inline const_reference at(const uint32_t& position) const{ return(this->__iterators[position]); }

private:
    template <class intrinsic_primitive>
    inline uint32_t GetNative(const buffer_type& buffer, const uint32_t position) const{
    	return(*reinterpret_cast<const intrinsic_primitive* const>(&buffer.buffer_[position*sizeof(intrinsic_primitive)]));
    }

private:
    size_type           n_entries;
    pointer             __iterators;
};

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_H_ */
