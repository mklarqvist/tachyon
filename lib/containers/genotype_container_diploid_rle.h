#ifndef CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_RLE_H_
#define CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_RLE_H_

#include "genotype_container_interface.h"

namespace tachyon{
namespace containers{

template <class T>
class GenotypeContainerDiploidRLE : public GenotypeContainerInterface{
private:
	typedef GenotypeContainerInterface    parent_type;
    typedef GenotypeContainerDiploidRLE   self_type;
    typedef T                             value_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;

public:
    GenotypeContainerDiploidRLE();
    GenotypeContainerDiploidRLE(const char* const data, const U32 n_entries, const meta_type& meta_entry);
    ~GenotypeContainerDiploidRLE();

    // Element access
    inline reference at(const size_type& position){ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline const_reference at(const size_type& position) const{ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline reference operator[](const size_type& position){ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline const_reference operator[](const size_type& position) const{ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline reference front(void){ return(*reinterpret_cast<const T* const>(&this->__data[0])); }
    inline const_reference front(void) const{ return(*reinterpret_cast<const T* const>(&this->__data[0])); }
    inline reference back(void){ return(*reinterpret_cast<const T* const>(&this->__data[(this->n_entries - 1) * sizeof(value_type)])); }
    inline const_reference back(void) const{ return(*reinterpret_cast<const T* const>(&this->__data[(this->n_entries - 1) * sizeof(value_type)])); }

    yon_gt* GetObjects(const uint32_t n_samples){
    	yon_gt* x = new yon_gt;
    	x->shift = this->__meta.IsAnyGTMissing()    ? 2 : 1;
		x->add   = this->__meta.IsGTMixedPhasing()  ? 1 : 0;
		x->global_phase = this->__meta.GetControllerPhase();
		x->data = this->__data;
		x->n_i = this->n_entries;
		x->method = 1;
		x->p = sizeof(T);
		x->m = 2;
		x->n_s = n_samples;
		x->n_allele = this->__meta.n_alleles;
		x->ppa = nullptr;

		assert(this->__meta.n_base_ploidy == 2);

		return(x);
    }

    yon_gt* GetObjects(yon_gt_ppa& ppa){
    	yon_gt* x = new yon_gt;
		x->shift = this->__meta.IsAnyGTMissing()    ? 2 : 1;
		x->add   = this->__meta.IsGTMixedPhasing()  ? 1 : 0;
		x->global_phase = this->__meta.GetControllerPhase();
		x->data = this->__data;
		x->m = 2;
		x->p = sizeof(T);
		x->n_i = this->n_entries;
		x->method = 1;
		x->n_s = ppa.n_samples;
		x->ppa = &ppa;
		x->n_allele = this->__meta.n_alleles;

		assert(this->__meta.n_base_ploidy == 2);

		return(x);
    }
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
GenotypeContainerDiploidRLE<T>::GenotypeContainerDiploidRLE(){

}

template <class T>
GenotypeContainerDiploidRLE<T>::GenotypeContainerDiploidRLE(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
	parent_type(data, n_entries, n_entries*sizeof(value_type), meta_entry)
{

}

template <class T>
GenotypeContainerDiploidRLE<T>::~GenotypeContainerDiploidRLE(){  }

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_RLE_H_ */
