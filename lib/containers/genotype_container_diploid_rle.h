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

template <class T>
U32 GenotypeContainerDiploidRLE<T>::getSum(void) const{
	U32 count = 0;
	const BYTE shift = this->__meta.IsAnyGTMissing()    ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing()  ? 1 : 0;

	for(U32 i = 0; i < this->n_entries; ++i)
		count += YON_GT_RLE_LENGTH(this->at(i), shift, add);

	return(count);
}

template <class T>
math::SquareMatrix<double>& GenotypeContainerDiploidRLE<T>::comparePairwise(square_matrix_type& square_matrix) const{
	// Has to be a SNV
	if(this->getMeta().IsBiallelicSNV() == false){
		//std::cerr << "skipping" << std::endl;
		return square_matrix;
	}

	const BYTE shift = this->__meta.IsAnyGTMissing()    ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing()  ? 1 : 0;

	U32 start_position = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		// self check
		const U32  ref_length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		const BYTE ref_alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		const BYTE ref_alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);

		// Cycle over implicit elements in object
		for(U32 start_sample = start_position; start_sample < start_position + ref_length; ++start_sample){
			for(U32 end_sample = start_sample + 1; end_sample < start_position + ref_length; ++end_sample){
				square_matrix(start_sample, end_sample) += 2;
				//++cum_length;
			}
		}
		//std::cerr << "(" << start_position << "," << start_position + ref_length << ")(" << start_position << "," << start_position + ref_length << ")\n";

		U32 internal_start = start_position + ref_length;

		// Compare to next object
		for(U32 j = i + 1; j < this->n_entries; ++j){
			const U32  length  = YON_GT_RLE_LENGTH(this->at(j), shift, add);
			const BYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(j), shift, add);
			const BYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(j), shift, add);
			//const U16 comp_genotype = (((this->at(j) & ((1 << shift) - 1) << add) >> add) << 8) | ((this->at(j) & ((1 << shift) - 1) << (add+shift)) >> (add+shift));
			const float score = this->comparatorSamplesDiploid(alleleA, ref_alleleA, alleleB, ref_alleleB);
			if(score == 0){
				internal_start += length;
				continue;
			}

			// Cycle over implicit elements in object
			//std::cerr << "(" << start_position << "," << start_position + ref_length << ")(" << internal_start << "," << internal_start + length << ")\n";
			for(U32 start_sample = start_position; start_sample < start_position + ref_length; ++start_sample){
				for(U32 end_sample = internal_start; end_sample < internal_start + length; ++end_sample){
					square_matrix(start_sample, end_sample) += score;
					//++cum_length;
				}
			}
			internal_start += length;
		}
		start_position += ref_length;
	}
	//std::cerr << start_position << std::endl;
	return(square_matrix);
}

/*
template <class T>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<T>::getLiteralObjects(void) const{
	std::vector<tachyon::core::GTObject> ret(this->n_entries);
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&ret[0]);
	for(U32 i = 0; i < this->n_entries; ++i)
		entries[i](this->at(i), this->__meta);

	return(ret);
}

template <class T>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<T>::getObjects(const U64& n_samples) const{
	std::vector<tachyon::core::GTObject> ret(n_samples);
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&ret[0]);

	const BYTE shift = this->__meta.IsAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing() ? 1 : 0;

	U32 cum_pos = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		const U32  length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		const BYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		const BYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);

		BYTE phasing = 0;
		if(add) phasing = this->at(i) & 1;
		else    phasing = this->__meta.GetControllerPhase();

		for(U32 j = 0; j < length; ++j, cum_pos++){
			entries[cum_pos].alleles = new core::GTObjectAllele[2];
			entries[cum_pos].alleles[0].allele  = alleleA;
			entries[cum_pos].alleles[1].allele  = alleleB;
			entries[cum_pos].alleles[0].phase = phasing;
			entries[cum_pos].alleles[1].phase = phasing;
			entries[cum_pos].n_objects = 1;
			entries[cum_pos].n_ploidy = 2;
		}
	}
	//assert(cum_pos == n_samples);
	return(ret);
}

template <class T>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<T>::getObjects(const U64& n_samples, const permutation_type& ppa_manager) const{
	std::vector<tachyon::core::GTObject> ret(n_samples);
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&ret[0]);

	const BYTE shift = this->__meta.IsAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing() ? 1 : 0;

	U32 cum_pos = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		const U32  length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		SBYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		SBYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);
		if(alleleA == 2) alleleA = -2;
		if(alleleB == 2) alleleB = -2;

		BYTE phasing = 0;
		if(add) phasing = this->at(i) & 1;
		else    phasing = this->__meta.GetControllerPhase();

		for(U32 j = 0; j < length; ++j, cum_pos++){
			entries[ppa_manager[cum_pos]].alleles = new core::GTObjectAllele[2];
			entries[ppa_manager[cum_pos]].alleles[0].allele  = alleleA;
			entries[ppa_manager[cum_pos]].alleles[1].allele  = alleleB;
			entries[ppa_manager[cum_pos]].alleles[0].phase = phasing;
			entries[ppa_manager[cum_pos]].alleles[1].phase = phasing;
			entries[ppa_manager[cum_pos]].n_objects = 1;
			entries[ppa_manager[cum_pos]].n_ploidy = 2;
		}
	}

	//assert(cum_pos == n_samples);
	return(ret);
}

template <class T>
void GenotypeContainerDiploidRLE<T>::getLiteralObjects(std::vector<tachyon::core::GTObject>& objects) const{
	if(objects.size() < this->size()) objects.resize(this->size());
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&objects[0]);

	const BYTE shift = this->__meta.IsAnyGTMissing() + 1;

	if(this->__meta.IsGTMixedPhasing()){
		for(U32 i = 0; i < this->size(); ++i)
			yon_gt_rle_obj<T>::Evaluate(entries[i], shift, this->at(i));

	} else {
		for(U32 i = 0; i < this->size(); ++i)
			yon_gt_rle_obj<T>::Evaluate(entries[i], shift, this->__meta.GetControllerPhase(), this->at(i));
	}
}

// Todo
template <class T>
void GenotypeContainerDiploidRLE<T>::getObjects(const U64& n_samples, std::vector<core::GTObject>& objects) const{
	if(objects.size() != n_samples){
		objects.resize(n_samples);
		for(U32 i = 0; i < n_samples; ++i)
			objects[i].alleles = new core::GTObjectAllele[2];
	}
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&objects[0]);

	const BYTE shift = this->__meta.IsAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing() ? 1 : 0;

	U32 cum_pos = 0;
	S32 alleleA, alleleB;
	U32 length;
	BYTE phasing = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);
		//alleleA -= core::YON_GT_RLE_CORRECTION[alleleA];
		//alleleB -= core::YON_GT_RLE_CORRECTION[alleleB];
		alleleA  = core::YON_GT_RLE_RECODE[alleleA];
		alleleB  = core::YON_GT_RLE_RECODE[alleleB];

		if(add) phasing = this->at(i) & 1;
		else    phasing = this->__meta.GetControllerPhase();

		for(U32 j = 0; j < length; ++j, cum_pos++){
			entries[cum_pos].alleles[0].allele  = alleleA;
			entries[cum_pos].alleles[1].allele  = alleleB;
			entries[cum_pos].alleles[0].phase = phasing;
			entries[cum_pos].alleles[1].phase = phasing;
			entries[cum_pos].n_objects = 1;
			entries[cum_pos].n_ploidy  = 2;
		}
	}
}

// Todo
template <class T>
void GenotypeContainerDiploidRLE<T>::getObjects(const U64& n_samples,
	                          std::vector<core::GTObject>& objects,
	                               const permutation_type& ppa_manager) const{
	if(objects.size() != n_samples){
		objects.resize(n_samples);
		for(U32 i = 0; i < n_samples; ++i)
			objects[i].alleles = new core::GTObjectAllele[2];
	}
	core::GTObjectDiploidRLE* entries = reinterpret_cast<core::GTObjectDiploidRLE*>(&objects[0]);

	const BYTE shift = this->__meta.IsAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta.IsGTMixedPhasing() ? 1 : 0;

	U32 cum_pos = 0;
	S32 alleleA, alleleB;
	U32 length;
	BYTE phasing = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);
		//alleleA -= core::YON_GT_RLE_CORRECTION[alleleA];
		//alleleB -= core::YON_GT_RLE_CORRECTION[alleleB];
		alleleA  = core::YON_GT_RLE_RECODE[alleleA];
		alleleB  = core::YON_GT_RLE_RECODE[alleleB];

		if(add) phasing = this->at(i) & 1;
		else    phasing = this->__meta.GetControllerPhase();

		for(U32 j = 0; j < length; ++j, cum_pos++){
			entries[ppa_manager[cum_pos]].alleles[0].allele  = alleleA;
			entries[ppa_manager[cum_pos]].alleles[1].allele  = alleleB;
			entries[ppa_manager[cum_pos]].alleles[0].phase = phasing;
			entries[ppa_manager[cum_pos]].alleles[1].phase = phasing;
			entries[ppa_manager[cum_pos]].n_objects = 1;
			entries[ppa_manager[cum_pos]].n_ploidy  = 2;
		}
	}
}
*/

template <class T>
GenotypeSummary& GenotypeContainerDiploidRLE<T>::updateSummary(gt_summary& gt_summary_object) const{
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
GenotypeSummary GenotypeContainerDiploidRLE<T>::getSummary(void) const{
	gt_summary summary;
	summary += *this;
	return(summary);
}

template <class T>
GenotypeSummary& GenotypeContainerDiploidRLE<T>::getSummary(gt_summary& gt_summary_object) const{
	gt_summary_object.clear();
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class T>
void GenotypeContainerDiploidRLE<T>::getTsTv(std::vector<ts_tv_object_type>& objects) const{
	if(this->size() == 0) return;
	if(this->getMeta().IsDiploid() == false) return;
	if(this->getMeta().alleles[0].size() != 1) return;
	assert(this->getMeta().GetNumberAlleles() == 2);

	BYTE references[10];
	switch(this->getMeta().alleles[0].allele[0]){
	case('A'): references[0] = constants::REF_ALT_A; break;
	case('T'): references[0] = constants::REF_ALT_T; break;
	case('G'): references[0] = constants::REF_ALT_G; break;
	case('C'): references[0] = constants::REF_ALT_C; break;
	default:   references[0] = constants::REF_ALT_MISSING; break;
	}

	if(references[0] == constants::REF_ALT_MISSING){
		std::cerr << "ref cannot be 0" << std::endl;
		return;
	}

	switch(this->getMeta().alleles[1].allele[0]){
	case('A'): references[1] = constants::REF_ALT_A; break;
	case('T'): references[1] = constants::REF_ALT_T; break;
	case('G'): references[1] = constants::REF_ALT_G; break;
	case('C'): references[1] = constants::REF_ALT_C; break;
	case('N'):
	default:   references[1] = constants::REF_ALT_MISSING; break;
	}

	references[2] = 4;
	references[3] = 4;
	references[4] = 4;
	references[5] = 4;
	references[6] = 4;
	references[7] = 4;
	references[8] = 4;
	references[9] = 4;

	const BYTE* const transition_map_target   = constants::TRANSITION_MAP[references[0]];
	const BYTE* const transversion_map_target = constants::TRANSVERSION_MAP[references[0]];
	const BYTE shift    = this->__meta.IsAnyGTMissing()   ? 2 : 1;
	const BYTE add      = this->__meta.IsGTMixedPhasing() ? 1 : 0;

	// Cycle over genotype objects
	U32 cum_position = 0;
	for(U32 i = 0; i < this->size(); ++i){
		const U32  length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		const BYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		const BYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);
		const BYTE targetA = references[alleleA];
		const BYTE targetB = references[alleleB];
		//const BYTE insertion_add = (references[alleleA] == constants::REF_ALT_INSERTION) + (references[alleleB] == constants::REF_ALT_INSERTION);
		const BYTE insertion_add = 0;

		for(U32 j = 0; j < length; ++j, cum_position++){
			objects[cum_position].n_insertions += insertion_add;
			++objects[cum_position].base_conversions[references[0]][targetA];
			++objects[cum_position].base_conversions[references[0]][targetB];
			objects[cum_position].n_transitions   += transition_map_target[targetA];
			objects[cum_position].n_transversions += transversion_map_target[targetA];
			objects[cum_position].n_transitions   += transition_map_target[targetB];
			objects[cum_position].n_transversions += transversion_map_target[targetB];
		}
	}
}

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_RLE_H_ */
