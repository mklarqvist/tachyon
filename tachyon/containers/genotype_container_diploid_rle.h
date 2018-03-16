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

    void operator()(const char* const data, const U32 n_entries, const meta_type& meta_entry){
		this->n_entries = n_entries;
		delete [] this->__data;

		const T* const re = reinterpret_cast<const T* const>(data);
		for(U32 i = 0; i < n_entries; ++i)
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

    // GT-specific
    U32 getSum(void) const;
    square_matrix_type& comparePairwise(square_matrix_type& square_matrix) const;
    std::vector<gt_object> getLiteralObjects(void) const;
    std::vector<gt_object> getObjects(const U64& n_samples) const;
    std::vector<gt_object> getObjects(const U64& n_samples, const permutation_type& ppa_manager) const;
    gt_summary& updateSummary(gt_summary& gt_summary_object) const;
    gt_summary getSummary(void) const;
    gt_summary& getSummary(gt_summary& gt_summary_object) const;
    void getTsTv(std::vector<ts_tv_object_type>& objects) const;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
GenotypeContainerDiploidRLE<return_type>::GenotypeContainerDiploidRLE(){

}

template <class return_type>
GenotypeContainerDiploidRLE<return_type>::GenotypeContainerDiploidRLE(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
	parent_type(data, n_entries, n_entries*sizeof(value_type), meta_entry)
{

}

template <class return_type>
GenotypeContainerDiploidRLE<return_type>::~GenotypeContainerDiploidRLE(){  }

template <class return_type>
U32 GenotypeContainerDiploidRLE<return_type>::getSum(void) const{
	U32 count = 0;
	const BYTE shift = this->__meta->isAnyGTMissing()    ? 2 : 1;
	const BYTE add   = this->__meta->isGTMixedPhasing()  ? 1 : 0;

	for(U32 i = 0; i < this->n_entries; ++i)
	count += YON_GT_RLE_LENGTH(this->at(i), shift, add);

	return(count);
}

template <class return_type>
math::SquareMatrix<double>& GenotypeContainerDiploidRLE<return_type>::comparePairwise(square_matrix_type& square_matrix) const{
	// Has to be a SNV
	if(this->getMeta().isSimpleSNV() == false){
		//std::cerr << "skipping" << std::endl;
		return square_matrix;
	}

	const BYTE shift = this->__meta->isAnyGTMissing()    ? 2 : 1;
	const BYTE add   = this->__meta->isGTMixedPhasing()  ? 1 : 0;

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

template <class return_type>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<return_type>::getLiteralObjects(void) const{
		std::vector<tachyon::core::GTObject> ret(this->n_entries);
		tachyon::core::GTObjectDiploidRLE* entries = reinterpret_cast<tachyon::core::GTObjectDiploidRLE*>(&ret[0]);
		for(U32 i = 0; i < this->n_entries; ++i)
			entries[i](this->at(i), *this->__meta);

		return(ret);
}

template <class return_type>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<return_type>::getObjects(const U64& n_samples) const{
	std::vector<tachyon::core::GTObject> ret(n_samples);
	tachyon::core::GTObjectDiploidRLE* entries = reinterpret_cast<tachyon::core::GTObjectDiploidRLE*>(&ret[0]);

	const BYTE shift = this->__meta->isAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta->isGTMixedPhasing() ? 1 : 0;

	U32 cum_pos = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		const U32  length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		const BYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		const BYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);

		BYTE phasing = 0;
		if(add) phasing = this->at(i) & 1;
		else    phasing = this->__meta->getControllerPhase();

		for(U32 j = 0; j < length; ++j, cum_pos++){
			entries[cum_pos].alleles = new std::pair<char,char>[2];
			entries[cum_pos].alleles[0].first  = alleleA;
			entries[cum_pos].alleles[1].first  = alleleB;
			entries[cum_pos].alleles[0].second = phasing;
			entries[cum_pos].alleles[1].second = phasing;
			entries[cum_pos].n_objects = 1;
			entries[cum_pos].n_alleles = 2;
		}
	}
	return(ret);
}

template <class return_type>
std::vector<tachyon::core::GTObject> GenotypeContainerDiploidRLE<return_type>::getObjects(const U64& n_samples, const permutation_type& ppa_manager) const{
	std::vector<gt_object> ret(this->getObjects(n_samples));
	std::vector<gt_object> ret_unpermuted(n_samples);

	for(U32 i = 0; i < n_samples; ++i)
		ret_unpermuted[i] = ret[ppa_manager[i]];

	//delete [] pointer;
	return(ret_unpermuted);
}

template <class return_type>
GenotypeSum& GenotypeContainerDiploidRLE<return_type>::updateSummary(gt_summary& gt_summary_object) const{
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class return_type>
GenotypeSum GenotypeContainerDiploidRLE<return_type>::getSummary(void) const{
	gt_summary summary;
	summary += *this;
	return(summary);
}

template <class return_type>
GenotypeSum& GenotypeContainerDiploidRLE<return_type>::getSummary(gt_summary& gt_summary_object) const{
	gt_summary_object += *this;
	return(gt_summary_object);
}

template <class return_type>
void GenotypeContainerDiploidRLE<return_type>::getTsTv(std::vector<ts_tv_object_type>& objects) const{
	std::cerr << "not yet" << std::endl;
	exit(1);

	if(this->size() == 0)
		return;

	// Has to be a SNV
	if(this->getMeta().isSimpleSNV() == false){
		//std::cerr << "skipping" << std::endl;
		return;
	}

	if(this->getMeta().isBiallelic() == false){
		return;
	}

	const BYTE shift = this->__meta->isAnyGTMissing()   ? 2 : 1;
	const BYTE add   = this->__meta->isGTMixedPhasing() ? 1 : 0;

	// If alleleA/B == ref then update self
	// If allele != ref then update ref->observed
	BYTE references[4];
	//references[0] = this->getMeta().getBiallelicAlleleLiteral(0);
	//references[1] = this->getMeta().getBiallelicAlleleLiteral(1);
	references[0] = 0;
	references[1] = 0;
	references[2] = 4; // Missing
	references[3] = 4; // EOV

	const BYTE* const transition_map_target   = constants::TRANSITION_MAP[references[0]];
	const BYTE* const transversion_map_target = constants::TRANSVERSION_MAP[references[0]];

	// Cycle over genotype objects
	U32 cum_position = 0;
	for(U32 i = 0; i < this->size(); ++i){
		const U32  length  = YON_GT_RLE_LENGTH(this->at(i), shift, add);
		const BYTE alleleA = YON_GT_RLE_ALLELE_A(this->at(i), shift, add);
		const BYTE alleleB = YON_GT_RLE_ALLELE_B(this->at(i), shift, add);

		if(alleleA == 2 || alleleB == 2){
			cum_position += length;
			continue;
		}

		for(U32 j = 0; j < length; ++j, cum_position++){
			++objects[cum_position].base_conversions[references[0]][references[alleleA]];
			++objects[cum_position].base_conversions[references[0]][references[alleleB]];
			objects[cum_position].n_transitions   += transition_map_target[references[alleleA]];
			objects[cum_position].n_transversions += transversion_map_target[references[alleleA]];
			objects[cum_position].n_transitions   += transition_map_target[references[alleleB]];
			objects[cum_position].n_transversions += transversion_map_target[references[alleleB]];
		}
	}
}

}
}



#endif /* CONTAINERS_GENOTYPE_CONTAINER_DIPLOID_RLE_H_ */
