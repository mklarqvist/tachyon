#ifndef CONTAINERS_GENOTYPE_CONTAINER_NPLOID_H_
#define CONTAINERS_GENOTYPE_CONTAINER_NPLOID_H_

#include "genotype_container_interface.h"

namespace tachyon{
namespace containers{

template <class T>
class GenotypeContainerNploid : public GenotypeContainerInterface{
private:
	typedef GenotypeContainerInterface    parent_type;
    typedef GenotypeContainerNploid       self_type;
    typedef T                             value_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;

public:
    GenotypeContainerNploid();
    GenotypeContainerNploid(const char* const data, const U32 n_entries, const meta_type& meta_entry);
    ~GenotypeContainerNploid();

    // Element access
    inline reference at(const size_type& position){ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline const_reference at(const size_type& position) const{ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline reference operator[](const size_type& position){ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline const_reference operator[](const size_type& position) const{ return(*reinterpret_cast<const T* const>(&this->__data[position * sizeof(value_type)])); }
    inline reference front(void){ return(*reinterpret_cast<const T* const>(&this->__data[0])); }
    inline const_reference front(void) const{ return(*reinterpret_cast<const T* const>(&this->__data[0])); }
    inline reference back(void){ return(*reinterpret_cast<const T* const>(&this->__data[(this->n_entries - 1) * sizeof(value_type)])); }
    inline const_reference back(void) const{ return(*reinterpret_cast<const T* const>(&this->__data[(this->n_entries - 1) * sizeof(value_type)])); }

    // GT-specific
    U32 getSum(void) const;
    square_matrix_type& comparePairwise(square_matrix_type& square_matrix) const;
    //std::vector<gt_object> getLiteralObjects(void) const;
    //std::vector<gt_object> getObjects(const U64& n_samples) const;
    //std::vector<gt_object> getObjects(const U64& n_samples, const permutation_type& ppa_manager) const;
    //void getObjects(const U64& n_samples, std::vector<gt_object>& objects) const;
	//void getObjects(const U64& n_samples, std::vector<gt_object>& objects, const permutation_type& ppa_manager) const;
	//void getLiteralObjects(std::vector<gt_object>& objects) const;

    gt_summary& updateSummary(gt_summary& gt_summary_object) const;
    gt_summary getSummary(void) const;
    gt_summary& getSummary(gt_summary& gt_summary_object) const;
    void getTsTv(std::vector<ts_tv_object_type>& objects) const;

    yon_gt* GetObjects(const uint32_t n_samples){
    	yon_gt* x = new yon_gt;
    	x->shift = 0;
    	x->add   = this->__meta.IsGTMixedPhasing();
		x->global_phase = this->__meta.GetControllerPhase();
		x->data = this->__data;
		x->n_i = this->n_entries;
		x->method = 4;
		x->m = this->__meta.n_base_ploidy;
		x->p = sizeof(T);
		x->n_s = n_samples;
		x->n_allele = this->__meta.n_alleles;
		x->ppa = nullptr;

		return(x);
    }

    yon_gt* GetObjects(yon_gt_ppa& ppa){
    	yon_gt* x = new yon_gt;
		x->shift = 0;
		x->add   = this->__meta.IsGTMixedPhasing();
		x->global_phase = this->__meta.GetControllerPhase();
		x->data = this->__data;
		x->m = this->__meta.n_base_ploidy;
		x->p = sizeof(T);
		x->n_i = this->n_entries;
		x->method = 4;
		x->n_s = ppa.n_samples;
		x->ppa = &ppa;
		x->n_allele = this->__meta.n_alleles;

		return(x);
    }
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
GenotypeContainerNploid<T>::GenotypeContainerNploid(){

}

template <class T>
GenotypeContainerNploid<T>::GenotypeContainerNploid(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
	parent_type(data, n_entries, n_entries*(sizeof(value_type) + meta_entry.n_base_ploidy*sizeof(uint8_t)), meta_entry)
{

}

template <class T>
GenotypeContainerNploid<T>::~GenotypeContainerNploid(){  }

template <class T>
U32 GenotypeContainerNploid<T>::getSum(void) const{
	U32 count = 0;
	const BYTE shift = this->__meta.IsAnyGTMissing()    ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing()  ? 1 : 0;

	for(U32 i = 0; i < this->n_entries; ++i)
		count += YON_GT_RLE_LENGTH(this->at(i), shift, add);

	return(count);
}

template <class T>
math::SquareMatrix<double>& GenotypeContainerNploid<T>::comparePairwise(square_matrix_type& square_matrix) const{
	return(square_matrix);
}


template <class T>
GenotypeSummary& GenotypeContainerNploid<T>::updateSummary(gt_summary& gt_summary_object) const{
	//gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
GenotypeSummary GenotypeContainerNploid<T>::getSummary(void) const{
	gt_summary summary;
	//summary += *this;
	return(summary);
}

template <class T>
GenotypeSummary& GenotypeContainerNploid<T>::getSummary(gt_summary& gt_summary_object) const{
	gt_summary_object.clear();
	//gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
void GenotypeContainerNploid<T>::getTsTv(std::vector<ts_tv_object_type>& objects) const{
	return;
}

}
}

#endif /* CONTAINERS_GENOTYPE_CONTAINER_NPLOID_H_ */
