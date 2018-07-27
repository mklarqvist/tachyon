#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "io/bcf/bcf_reader.h"
#include "permutation_manager.h"
#include "core/genotype_summary.h"

#include "io/htslib_integration.h"

namespace tachyon {
namespace algorithm {

struct yon_radix_gt {
	yon_radix_gt(void) : n_id(0), n_allocated(0), n_alleles(0), alleles(nullptr){}
	yon_radix_gt(const uint64_t n_id, const uint8_t n_allocate) :
		n_id(n_id),
		n_allocated(n_allocate),
		n_alleles(0),
		alleles(new int32_t[n_allocate])
	{

	}
	yon_radix_gt(const yon_radix_gt& other) :
		n_id(other.n_id),
		n_allocated(other.n_allocated),
		n_alleles(other.n_alleles),
		alleles(new int32_t[other.n_allocated])
	{
		memcpy(this->alleles, other.alleles, sizeof(int32_t)*other.n_alleles);
	}

	yon_radix_gt& operator=(const yon_radix_gt& other){
		delete [] this->alleles;
		this->n_id        = other.n_id;
		this->n_allocated = other.n_allocated;
		this->n_alleles   = other.n_alleles;
		this->alleles     = new int32_t[other.n_allocated];
		memcpy(this->alleles, other.alleles, sizeof(int32_t)*other.n_alleles);
		return(*this);
	}

	~yon_radix_gt(){ delete [] this->alleles; }

	inline int32_t& operator[](const int32_t value){ return(this->alleles[value]); }
	inline const int32_t& operator[](const int32_t value) const{ return(this->alleles[value]); }

	void clear(void){
		this->n_id = 0;
		memset(this->alleles, 0, sizeof(int32_t)*this->n_alleles);
		this->n_alleles = 0;
	}

	void Allocate(const uint64_t n_id, const U32 n_allocate){
		delete [] this->alleles;
		this->n_id        = n_id;
		this->n_allocated = n_allocate;
		this->n_alleles   = 0;
		this->alleles     = new int32_t[this->n_allocated];
	}

	bool operator<(const yon_radix_gt& other) const {
		if(this->n_alleles < other.n_alleles) return true;
		if(other.n_alleles < this->n_alleles) return false;

		for(U32 i = 0; i < this->n_alleles; ++i){
			if(this->alleles[i] < other.alleles[i]){
				return true;
			}
		}
		return false;
	}

	inline bool operator!=(const yon_radix_gt& other) const{ return(!(*this == other)); }

	bool operator==(const yon_radix_gt& other) const{
		// Do not compare identifier.
		if(this->n_alleles != other.n_alleles)
			return false;

		for(U32 i = 0; i < this->n_alleles; ++i){
			if(this->alleles[i] != other.alleles[i])
				return false;
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const yon_radix_gt& gt){
		stream << gt.n_id << ": ";
		if(gt.n_alleles){
			stream << (uint32_t)gt.alleles[0];
			for(U32 i = 1; i < gt.n_alleles; ++i)
				stream << "," << (uint32_t)gt.alleles[i];
		}
		return stream;
	}

	// todo: sort
	uint64_t n_id;
	uint8_t  n_allocated;
	uint8_t  n_alleles;
	int32_t* alleles;
};

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * bi-allelic diploid.
 */
class RadixSortGT {
	typedef RadixSortGT        self_type;
	typedef bcf::BCFReader     bcf_reader_type;
	typedef bcf::BCFEntry      bcf_entry_type;
	typedef PermutationManager manager_type;

	typedef containers::VcfContainer vcf_container_type;

public:
	RadixSortGT();
	RadixSortGT(const U64 n_samples);
	~RadixSortGT();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	void reset(void);
	void SetSamples(const U64 n_samples);

	bool Build(const vcf_container_type& vcf_container);

	// Construct given a reader with a block
	// of BCF entries loaded in it
	bool Build(const bcf_reader_type& reader);
	bool Update(const bcf_entry_type& entry);

	inline const U64& GetNumberSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }

	/**<
	* Static function that calculates the 64-bit hash value for the target
	* FORMAT/FILTER/INFO vector of id fields. The id fields must be of type
	* int (S32). Example of using this function:
	*
	* const U64 hash_value = VariantImporter::HashIdentifiers(id_vector);
	*
	* @param id_vector Input vector of FORMAT/FILTER/INFO identifiers.
	* @return          Returns a 64-bit hash value.
	*/
	template <class T>
	static U64 HashGenotypes(const T* genotype_data, const U32 l_data) {
		XXH64_state_t* const state = XXH64_createState();
		if (state==NULL) exit(1);

		XXH_errorcode const resetResult = XXH64_reset(state, BCF_HASH_SEED);
		if (resetResult == XXH_ERROR) exit(1);

		for(U32 i = 0; i < l_data; ++i) {
			XXH_errorcode const addResult = XXH64_update(state, (const void*)genotype_data, sizeof(T));
			if (addResult == XXH_ERROR) exit(1);
		}

		U64 hash = XXH64_digest(state);
		XXH64_freeState(state);

		return hash;
	}

public:
	U64           n_samples; // total number of entries in file
	U32           position;  // number of entries parsed
	U32           p_i[9];    // number of entries in bin i
	BYTE*         GT_array;  // packed genotype array
	U32**         bins;      // bin i
	manager_type* manager;   // permutation manager
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
