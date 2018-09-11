#ifndef CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_
#define CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_

#include "genotype_container_interface.h"

namespace tachyon{
namespace containers{

template <class T>
class GenotypeContainerDiploidBCF : public GenotypeContainerInterface{
public:
	typedef GenotypeContainerInterface    parent_type;
    typedef GenotypeContainerDiploidBCF   self_type;
    typedef T                             value_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;

public:
    GenotypeContainerDiploidBCF();
    GenotypeContainerDiploidBCF(const char* const data,
		const size_type& n_entries,
		const uint16_t n_als,
		const uint8_t  n_ploidy,
		const uint16_t controller) :
			parent_type(data,n_entries,n_entries*sizeof(value_type),n_als,n_ploidy,controller)
	{

	}
    ~GenotypeContainerDiploidBCF();

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
		return(nullptr);
	}

	yon_gt* GetObjects(yon_gt_ppa& ppa){
		return(nullptr);
	}
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
GenotypeContainerDiploidBCF<T>::GenotypeContainerDiploidBCF(){

}

template <class T>
GenotypeContainerDiploidBCF<T>::~GenotypeContainerDiploidBCF(){  }

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_ */
