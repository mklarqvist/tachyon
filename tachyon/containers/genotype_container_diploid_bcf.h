#ifndef CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_
#define CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_

#include "genotype_container_interface.h"

namespace tachyon{
namespace containers{

template <class T>
class GenotypeContainerDiploidBCF : public GenotypeContainerInterface{
private:
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
    GenotypeContainerDiploidBCF(const char* const data, const U32 n_entries, const meta_type& meta_entry);
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

    // GT-specific
    U32 getSum(void) const;
    square_matrix_type& comparePairwise(square_matrix_type& square_matrix) const;
    std::vector<gt_object> getLiteralObjects(void) const;
    std::vector<gt_object> getObjects(const U64& n_samples) const;
    std::vector<gt_object> getObjects(const U64& n_samples, const permutation_type& ppa_manager) const;
    void getLiteralObjects(std::vector<gt_object>& objects) const;
	void getObjects(std::vector<gt_object>& objects, const U64& n_samples) const;
	void getObjects(std::vector<gt_object>& objects, const U64& n_samples, const permutation_type& ppa_manager) const;
    gt_summary& updateSummary(gt_summary& gt_summary_object) const;
    gt_summary getSummary(void) const;
    gt_summary& getSummary(gt_summary& gt_summary_object) const;
    void getTsTv(std::vector<ts_tv_object_type>& objects) const;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
GenotypeContainerDiploidBCF<T>::GenotypeContainerDiploidBCF(){

}

template <class T>
GenotypeContainerDiploidBCF<T>::GenotypeContainerDiploidBCF(const char* const  data,
                                                                    const U32  n_entries,
                                                              const meta_type& meta_entry) :
	parent_type(data, n_entries, n_entries*sizeof(value_type), meta_entry)
{

}

template <class T>
GenotypeContainerDiploidBCF<T>::~GenotypeContainerDiploidBCF(){  }

template <class T>
U32 GenotypeContainerDiploidBCF<T>::getSum(void) const{
	return(this->size());
}

template <class T>
math::SquareMatrix<double>& GenotypeContainerDiploidBCF<T>::comparePairwise(square_matrix_type& square_matrix) const{
	std::cerr << "not implemented" << std::endl;
	exit(1);
	return(square_matrix);
}

template <class T>
std::vector<core::GTObject> GenotypeContainerDiploidBCF<T>::getLiteralObjects(void) const{
	std::vector<core::GTObject> ret(this->n_entries);
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&ret[0]);
	for(U32 i = 0; i < this->n_entries; ++i)
		entries[i](this->at(i), this->__meta);

	return(ret);
}

template <class T>
std::vector<core::GTObject> GenotypeContainerDiploidBCF<T>::getObjects(const U64& n_samples) const{
	std::vector<core::GTObject> ret(n_samples);
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&ret[0]);

	const BYTE shift    = (sizeof(T)*8 - 1) / 2;

	for(U32 i = 0; i < this->n_entries; ++i){
		entries[i].alleles = new std::pair<char,char>[2];
		entries[i].alleles[0].first  = YON_GT_DIPLOID_BCF_A(this->at(i), shift);
		entries[i].alleles[1].first  = YON_GT_DIPLOID_BCF_B(this->at(i), shift);
		entries[i].alleles[0].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[i].alleles[1].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[i].n_objects = 1;
		entries[i].n_alleles = 2;
	}
	//assert(cum_pos == n_samples);
	return(ret);
}

template <class T>
std::vector<core::GTObject> GenotypeContainerDiploidBCF<T>::getObjects(const U64& n_samples, const permutation_type& ppa_manager) const{
	std::vector<core::GTObject> ret(n_samples);
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&ret[0]);

	const BYTE shift = (sizeof(T)*8 - 1) / 2;

	U32 cum_pos = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		entries[ppa_manager[cum_pos]].alleles = new std::pair<char,char>[2];
		entries[ppa_manager[cum_pos]].alleles[0].first  = YON_GT_DIPLOID_BCF_A(this->at(i), shift);
		entries[ppa_manager[cum_pos]].alleles[1].first  = YON_GT_DIPLOID_BCF_B(this->at(i), shift);
		entries[ppa_manager[cum_pos]].alleles[0].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[ppa_manager[cum_pos]].alleles[1].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[ppa_manager[cum_pos]].n_objects = 1;
		entries[ppa_manager[cum_pos]].n_alleles = 2;
	}

	//assert(cum_pos == n_samples);
	return(ret);
}

template <class T>
void GenotypeContainerDiploidBCF<T>::getLiteralObjects(std::vector<core::GTObject>& objects) const{
	if(objects.size() < this->size()) objects.resize(this->size());
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&objects[0]);
	for(U32 i = 0; i < this->n_entries; ++i)
		entries[i](this->at(i), this->__meta);
}

template <class T>
void GenotypeContainerDiploidBCF<T>::getObjects(std::vector<core::GTObject>& objects, const U64& n_samples) const{
	if(objects.size() < n_samples) objects.resize(n_samples);
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&objects[0]);

	const BYTE shift    = (sizeof(T)*8 - 1) / 2;

	for(U32 i = 0; i < this->n_entries; ++i){
		delete [] entries[i].alleles;
		entries[i].alleles = new std::pair<char,char>[2];
		entries[i].alleles[0].first  = YON_GT_DIPLOID_BCF_A(this->at(i), shift);
		entries[i].alleles[1].first  = YON_GT_DIPLOID_BCF_B(this->at(i), shift);
		entries[i].alleles[0].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[i].alleles[1].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[i].n_objects = 1;
		entries[i].n_alleles = 2;
	}
}

template <class T>
void GenotypeContainerDiploidBCF<T>::getObjects(std::vector<core::GTObject>& objects, const U64& n_samples, const permutation_type& ppa_manager) const{
	if(objects.size() < n_samples) objects.resize(n_samples);
	core::GTObjectDiploidBCF* entries = reinterpret_cast<core::GTObjectDiploidBCF*>(&objects[0]);

	const BYTE shift = (sizeof(T)*8 - 1) / 2;

	U32 cum_pos = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		delete [] entries[ppa_manager[cum_pos]].alleles;
		entries[ppa_manager[cum_pos]].alleles = new std::pair<char,char>[2];
		entries[ppa_manager[cum_pos]].alleles[0].first  = YON_GT_DIPLOID_BCF_A(this->at(i), shift);
		entries[ppa_manager[cum_pos]].alleles[1].first  = YON_GT_DIPLOID_BCF_B(this->at(i), shift);
		entries[ppa_manager[cum_pos]].alleles[0].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[ppa_manager[cum_pos]].alleles[1].second = YON_GT_DIPLOID_BCF_PHASE(this->at(i));
		entries[ppa_manager[cum_pos]].n_objects = 1;
		entries[ppa_manager[cum_pos]].n_alleles = 2;
	}
}

template <class T>
GenotypeSummary& GenotypeContainerDiploidBCF<T>::updateSummary(gt_summary& gt_summary_object) const{
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
GenotypeSummary GenotypeContainerDiploidBCF<T>::getSummary(void) const{
	gt_summary summary;
	summary += *this;
	return(summary);
}

template <class T>
GenotypeSummary& GenotypeContainerDiploidBCF<T>::getSummary(gt_summary& gt_summary_object) const{
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
void GenotypeContainerDiploidBCF<T>::getTsTv(std::vector<ts_tv_object_type>& objects) const{
	if(this->size() == 0)
		return;

	// Has to be a SNV and biallelic
	if(this->getMeta().isBiallelicSNV() == false) return;

	const BYTE shift = (sizeof(T)*8 - 1) / 2;

	// If alleleA/B == ref then update self
	// If allele != ref then update ref->observed
	BYTE references[4];
	//references[0] = this->getMeta().getBiallelicAlleleLiteral(0);
	//references[1] = this->getMeta().getBiallelicAlleleLiteral(1);
	references[0] = 0;
	references[1] = 0;
	references[2] = 4; // Missing
	references[3] = 4; // EOV: is never available in this encodin

	const BYTE* const transition_map_target   = constants::TRANSITION_MAP[references[0]];
	const BYTE* const transversion_map_target = constants::TRANSVERSION_MAP[references[0]];

	// Cycle over genotype objects
	for(U32 i = 0; i < this->size(); ++i){
		const BYTE alleleA = YON_GT_DIPLOID_BCF_A(this->at(i), shift);
		const BYTE alleleB = YON_GT_DIPLOID_BCF_B(this->at(i), shift);

		++objects[i].base_conversions[references[0]][references[alleleA]];
		++objects[i].base_conversions[references[0]][references[alleleB]];
		objects[i].n_transitions   += transition_map_target[references[alleleA]];
		objects[i].n_transversions += transversion_map_target[references[alleleA]];
		objects[i].n_transitions   += transition_map_target[references[alleleB]];
		objects[i].n_transversions += transversion_map_target[references[alleleB]];
	}
}

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_BCF_H_ */
