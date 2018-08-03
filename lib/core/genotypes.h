#ifndef CORE_GENOTYPES_H_
#define CORE_GENOTYPES_H_

#include <cstring>

#include "io/basic_buffer.h"
#include "third_party/intervalTree.h"

namespace tachyon{

// Basic structure that maintains the permutation
// order of the samples in relation to the global header.
// This object is required if you want to use individual
// genotypes in the ORIGINAL order. If this is not required
// in your use-case then this structure has no value.
struct yon_gt_ppa {
	yon_gt_ppa(void) : n_samples(0), ordering(nullptr){}
	yon_gt_ppa(const uint32_t n_samples) : n_samples(n_samples), ordering(new uint32_t[n_samples]){ this->reset(); }
	~yon_gt_ppa(void){ delete [] this->ordering; }

	uint32_t& operator[](const uint32_t& position){ return(this->ordering[position]); }
	const uint32_t& operator[](const uint32_t& position) const{ return(this->ordering[position]); }
	uint32_t& at(const uint32_t& position){ return(this->ordering[position]); }
	const uint32_t& at(const uint32_t& position) const{ return(this->ordering[position]); }

	void Allocate(const uint32_t n_samples){
		delete [] this->ordering;
		this->n_samples = n_samples;
		this->ordering = new uint32_t[n_samples];
		this->reset();
	}

	void reset(void){
		for(U32 i = 0; i < this->n_samples; ++i)
			this->ordering[i] = i;
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, yon_gt_ppa& ppa){
		io::DeserializePrimitive(ppa.n_samples, buffer);
		ppa.ordering = new uint32_t[ppa.n_samples];
		for(U32 i = 0; i < ppa.n_samples; ++i)
			io::DeserializePrimitive(ppa.ordering[i], buffer);

		return(buffer);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const yon_gt_ppa& ppa){
		io::SerializePrimitive(ppa.n_samples, buffer);
		for(U32 i = 0; i < ppa.n_samples; ++i)
			io::SerializePrimitive(ppa.ordering[i], buffer);

		return(buffer);
	}

	uint32_t  n_samples;
	uint32_t* ordering;
};

struct yon_radix_gt {
	yon_radix_gt() : n_ploidy(0), n_allocated(4), id(0), alleles(new uint16_t[this->n_allocated]){
		memset(this->alleles, 0, sizeof(uint16_t)*this->n_allocated);
	}
	~yon_radix_gt(){ delete [] this->alleles; }
	yon_radix_gt(const yon_radix_gt& other) : n_ploidy(other.n_ploidy), n_allocated(other.n_allocated), id(other.id), alleles(new uint16_t[this->n_allocated])
	{
		memcpy(this->alleles, other.alleles, sizeof(uint16_t)*this->n_allocated);
	}
	yon_radix_gt(yon_radix_gt&& other) : n_ploidy(other.n_ploidy), n_allocated(other.n_allocated), id(other.id), alleles(other.alleles)
	{
		other.alleles = nullptr;
	}

	yon_radix_gt& operator=(const yon_radix_gt& other) // copy assignment
	{
		this->id = other.id;
		this->n_ploidy = other.n_ploidy;
		this->n_allocated = other.n_allocated;
		delete [] this->alleles;
		this->alleles = new uint16_t[this->n_allocated];
		memcpy(this->alleles, other.alleles, sizeof(uint16_t)*this->n_allocated);
		return *this;
	}
	yon_radix_gt& operator=(yon_radix_gt&& other) // move assignment
	{
		if(this!=&other) // prevent self-move
		{
			this->id = other.id;
			this->n_ploidy = other.n_ploidy;
			this->n_allocated = other.n_allocated;
			delete [] this->alleles;
			this->alleles = other.alleles;
			other.alleles = nullptr;
		}
		return *this;
	}

	bool operator<(const yon_radix_gt& other) const{
		// Do not compare incremental sample identification
		// numbers as that is not the desired outcome of
		// the sort.
		if(this->n_ploidy < other.n_ploidy) return true;
		if(other.n_ploidy < this->n_ploidy) return false;

		for(U32 i = 0; i < this->n_ploidy; ++i){
			if(this->alleles[i] < other.alleles[i])
				return true;
		}
		return false;
	}

	bool operator==(const yon_radix_gt& other) const{
		// Do not compare incremental sample identification
		// numbers as that is not the desired outcome of
		// the comparison.
		if(this->n_ploidy != other.n_ploidy)
			return false;

		for(U32 i = 0; i < this->n_ploidy; ++i){
			if(this->alleles[i] != other.alleles[i])
				return false;
		}
		return true;
	}

	inline bool operator!=(const yon_radix_gt& other) const{ return(!(*this == other)); }

	friend std::ostream& operator<<(std::ostream& stream, const yon_radix_gt& genotype){
		stream << genotype.id << ":";
		if(genotype.n_ploidy){
			stream << genotype.alleles[0];
			for(U32 i = 1; i < genotype.n_ploidy; ++i){
				stream << "," << genotype.alleles[i];
			}
		}
		return(stream);
	}

	U64 GetPackedInteger(const uint8_t& shift_size = 8) const{
		U64 packed = 0;
		for(U32 i = 0; i < this->n_ploidy; ++i){
			packed <<= shift_size;
			assert(((this->alleles[i] << shift_size) >> shift_size) == this->alleles[i]);
			packed |= (this->alleles[i] & ((1 << shift_size)) - 1);
		}
		return packed;
	}

	void resize(const uint8_t new_ploidy){
		uint16_t* temp = new uint16_t[new_ploidy];
		memcpy(temp, this->alleles, this->n_allocated * sizeof(uint16_t));
		delete [] this->alleles;
		this->alleles = temp;
		this->n_allocated = new_ploidy;
	}

	uint8_t   n_ploidy;
	uint8_t   n_allocated;
	uint64_t  id;
	uint16_t* alleles;
};

struct yon_gt_rcd {
	uint32_t length : 31, phase : 1;
	//uint8_t alleles; // Determined from base ploidy
	uint8_t* allele;
};

struct yon_gt {
    uint8_t  add : 7,
             global_phase : 1;
    uint8_t  shift;
    uint8_t  p, m, method; // bytes per entry, base ploidy, base method
    uint32_t n_s, n_i;     // number samples, number of entries
    yon_gt_ppa* ppa; // pointer to ppa
    uint8_t* data; // pointer to data
    // Todo: Lazy evaluate (d)
    // Todo: if ppa then map from sample -> internal record
    uint8_t* d_bcf; // lazy evaluated as Bcf entries (length = base_ploidy * n_samples * sizeof(uint8_t))
    uint8_t* d_ppa; // lazy evaluated from ppa to internal offset (length = n_samples)
    yon_gt_rcd* rcds; // lazy interpreted internal records
    // Todo: interval tree of cumulative runs -> ppa find overlap
    // e.g. itree[45] -> intersect overlap
    algorithm::IntervalTree<uint32_t, uint32_t>* itree; // interval tree for consecutive ranges

    bool EvaluateBcf();
    bool EvaluatePpa();
    bool EvaluateRecords();

    template <class T>
    inline T& GetPrimitive(const uint32_t sample){ return(*reinterpret_cast<T*>(&this->data[sample])); }

    template <class T>
	inline T& GetPrimitivePPA(const uint32_t sample){ return(this->ppa[*reinterpret_cast<T*>(&this->data[sample])]); }
};

}



#endif /* CORE_GENOTYPES_H_ */
