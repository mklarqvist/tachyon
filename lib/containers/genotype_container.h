#ifndef CONTAINERS_GENOTYPE_CONTAINER_H_
#define CONTAINERS_GENOTYPE_CONTAINER_H_

#include "genotype_container_diploid_bcf.h"
#include "genotype_container_diploid_rle.h"
#include "genotype_container_diploid_simple.h"
#include "genotype_container_nploid.h"
#include "meta_container.h"
#include "variant_block.h"

namespace tachyon{
namespace containers{

// forward declaration for getters
template <class T> class GenotypeContainerDiploidRLE;
template <class T> class GenotypeContainerDiploidSimple;


class GenotypeContainer{
private:
    typedef GenotypeContainer          self_type;
    typedef GenotypeContainerInterface value_type;
    typedef value_type&                reference;
    typedef const value_type&          const_reference;
    typedef value_type*                pointer;
    typedef const value_type*          const_pointer;
    typedef std::ptrdiff_t             difference_type;
    typedef std::size_t                size_type;
    typedef MetaContainer              meta_container_type;
    typedef tachyon::core::MetaEntry   meta_type;
    typedef io::BasicBuffer            buffer_type;
    typedef containers::GenotypeSummary      gt_summary_type;
    typedef VariantBlock               block_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    GenotypeContainer(const block_type& block, const MetaContainer& meta);
    ~GenotypeContainer();

    // Capacity
	inline bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Element access
	inline pointer         data(void){ return(this->__iterators); }
	inline const_pointer   data(void) const{ return(this->__iterators); }
	inline reference       operator[](const U32& position){ return(this->__iterators[position]); }
	inline const_reference operator[](const U32& position) const{ return(this->__iterators[position]); }
	inline reference       at(const U32& position){ return(this->__iterators[position]); }
	inline const_reference at(const U32& position) const{ return(this->__iterators[position]); }

	// Advanced getters
	inline GenotypeContainerDiploidSimple<BYTE>* getDiploidSimpleByte(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidSimple<BYTE>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidSimple<U16>*  getDiploidSimpleU16(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidSimple<U16>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidSimple<U32>*  getDiploidSimpleU32(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidSimple<U32>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidSimple<U64>*  getDiploidSimpleU64(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidSimple<U64>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidRLE<BYTE>*    getDiploidRLEByte(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidRLE<BYTE>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidRLE<U16>*     getDiploidRLEU16(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidRLE<U16>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidRLE<U32>*     getDiploidRLEU32(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidRLE<U32>*>(&this->__iterators[position])); }
	inline GenotypeContainerDiploidRLE<U64>*     getDiploidRLEU64(const U32 position){ return(reinterpret_cast<GenotypeContainerDiploidRLE<U64>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidSimple<BYTE>* getDiploidSimpleByte(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidSimple<BYTE>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidSimple<U16>*  getDiploidSimpleU16(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidSimple<U16>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidSimple<U32>*  getDiploidSimpleU32(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidSimple<U32>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidSimple<U64>*  getDiploidSimpleU64(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidSimple<U64>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidRLE<BYTE>*    getDiploidRLEByte(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidRLE<BYTE>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidRLE<U16>*     getDiploidRLEU16(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidRLE<U16>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidRLE<U32>*     getDiploidRLEU32(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidRLE<U32>*>(&this->__iterators[position])); }
	inline const GenotypeContainerDiploidRLE<U64>*     getDiploidRLEU64(const U32 position) const{ return(reinterpret_cast<GenotypeContainerDiploidRLE<U64>*>(&this->__iterators[position])); }

private:
    template <class intrinsic_primitive> inline U32 getNative(const buffer_type& buffer, const U32 position) const{
    	return(*reinterpret_cast<const intrinsic_primitive* const>(&buffer.buffer[position*sizeof(intrinsic_primitive)]));
    }

private:
    size_type           n_entries;
    meta_container_type __meta_container;
    pointer             __iterators;
};

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_H_ */
