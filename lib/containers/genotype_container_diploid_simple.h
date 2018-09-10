#ifndef CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_SIMPLE_H_
#define CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_SIMPLE_H_

#include "genotype_container_interface.h"

namespace tachyon{
namespace containers{

template <class T>
class GenotypeContainerDiploidSimple : public GenotypeContainerInterface{
public:
	typedef GenotypeContainerInterface       parent_type;
    typedef GenotypeContainerDiploidSimple   self_type;
    typedef T                                value_type;
    typedef value_type&                      reference;
    typedef const value_type&                const_reference;
    typedef value_type*                      pointer;
    typedef const value_type*                const_pointer;
    typedef std::ptrdiff_t                   difference_type;
    typedef std::size_t                      size_type;

public:
    GenotypeContainerDiploidSimple();
    GenotypeContainerDiploidSimple(const char* const data, const uint32_t n_entries, const meta_type& meta_entry);
    GenotypeContainerDiploidSimple(const char* const data,
		const size_type& n_entries,
		const uint16_t n_als,
		const uint8_t  n_ploidy,
		const uint16_t controller) :
			parent_type(data,n_entries,n_entries*sizeof(value_type),n_als,n_ploidy,controller)
	{

	}
    ~GenotypeContainerDiploidSimple();

    void operator()(const char* const data, const uint32_t n_entries, const meta_type& meta_entry){
		this->n_entries = n_entries;
		delete [] this->__data;

		const T* const re = reinterpret_cast<const T* const>(data);
		for(uint32_t i = 0; i < n_entries; ++i)
			this->__data[i] = re[i];
    }

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
		x->n_allele = this->n_alleles;
		x->shift = ceil(log2(x->n_allele + 2 + 1));
		x->add   =  this->gt_mixed_phase ? 1 : 0;
		x->global_phase = this->gt_global_phase;
		x->data = this->__data;
		x->n_i = this->n_entries;
		x->method = 2;
		x->p = sizeof(T);
		x->m = 2;
		x->n_s = n_samples;
		x->ppa = nullptr;

		assert(n_base_ploidy == 2);

		return(x);
	}

	yon_gt* GetObjects(yon_gt_ppa& ppa){
		yon_gt* x = new yon_gt;
		x->n_allele = this->n_alleles;
		x->shift = ceil(log2(x->n_allele + 2 + 1));
		x->add   =  this->gt_mixed_phase ? 1 : 0;
		x->global_phase = this->gt_global_phase;
		x->data = this->__data;
		x->m = 2;
		x->p = sizeof(T);
		x->n_i = this->n_entries;
		x->method = 2;
		x->n_s = ppa.n_s;
		x->ppa = &ppa;

		assert(this->n_base_ploidy == 2);

		return(x);
	}
};



// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
GenotypeContainerDiploidSimple<return_type>::GenotypeContainerDiploidSimple(){

}

template <class return_type>
GenotypeContainerDiploidSimple<return_type>::GenotypeContainerDiploidSimple(const char* const data, const uint32_t n_entries, const meta_type& meta_entry) :
	parent_type(data, n_entries, n_entries*sizeof(value_type), meta_entry)
{

}

template <class return_type>
GenotypeContainerDiploidSimple<return_type>::~GenotypeContainerDiploidSimple(){  }

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_SIMPLE_H_ */
