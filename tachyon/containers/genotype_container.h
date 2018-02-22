#ifndef CONTAINERS_GENOTYPE_CONTAINER_H_
#define CONTAINERS_GENOTYPE_CONTAINER_H_

#include "datablock.h"
#include "genotype_container_diploid_rle.h"
#include "genotype_container_diploid_simple.h"
#include "meta_container.h"

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
    typedef containers::GenotypeSum    gt_summary_type;

    // Function pointers
	typedef const U32 (self_type::*getNativeFuncDef)(const buffer_type& buffer, const U32 position) const;

public:
    GenotypeContainer(const DataBlock& block);
    ~GenotypeContainer();

    class iterator{
        private:
    		typedef iterator self_type;
    		typedef std::forward_iterator_tag iterator_category;

        public:
    		iterator(pointer ptr) : ptr_(ptr) { }
    		void operator++() { ptr_++; }
    		void operator++(int junk) { ptr_++; }
    		reference operator*() const{ return *ptr_; }
    		pointer operator->() const{ return ptr_; }
    		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
    		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
    	private:
    		pointer ptr_;
    	};

        class const_iterator{
    	private:
    		typedef const_iterator self_type;
    		typedef std::forward_iterator_tag iterator_category;

    	public:
    		const_iterator(pointer ptr) : ptr_(ptr) { }
    		void operator++() { ptr_++; }
    		void operator++(int junk) { ptr_++; }
    		const_reference operator*() const{ return *ptr_; }
    		const_pointer operator->() const{ return ptr_; }
    		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
    		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
    	private:
    		pointer ptr_;
    	};

    // Capacity
	inline const bool       empty(void) const{ return(this->n_entries == 0); }
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

	// GT summary
	gt_summary_type calculateSummary(void) const{
		gt_summary_type summary;
		for(U32 i = 0; i < this->size(); ++i){
			this->at(i).updateSummary(summary);
		}
		return(summary);
	}


private:
    template <class intrinsic_primitive> inline const U32 getNative(const buffer_type& buffer, const U32 position) const{
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
