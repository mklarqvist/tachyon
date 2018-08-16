#ifndef CORE_GENOTYPES_H_
#define CORE_GENOTYPES_H_

#include <cstring>

#include "htslib/vcf.h"

#include "io/basic_buffer.h"
#include "containers/components/generic_iterator.h"
#include "third_party/intervalTree.h"
#include "math/fisher_math.h"
#include "utility/support_vcf.h"

namespace tachyon{

#define YON_GT_RLE_ALLELE_A(PRIMITITVE, SHIFT, ADD)  (((PRIMITITVE) & ((1 << (SHIFT)) - 1) << (ADD)) >> (ADD));
#define YON_GT_RLE_ALLELE_B(PRIMITIVE, SHIFT, ADD)   (((PRIMITIVE) & ((1 << (SHIFT)) - 1) << ((ADD)+(SHIFT))) >> ((ADD)+(SHIFT)));
#define YON_GT_RLE_LENGTH(PRIMITIVE, SHIFT, ADD)     ((PRIMITIVE) >> (2*(SHIFT) + (ADD)))
#define YON_GT_DIPLOID_ALLELE_LOOKUP(A,B,shift,mask) (((A) & (mask)) << (shift)) | ((B) & (mask))
#define YON_GT_DIPLOID_BCF_A(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> ((SHIFT) + 1)) & (((U64)1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_B(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> 1) & (((U64)1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_PHASE(PRIMITIVE)          ((PRIMITIVE) & 1)
#define YON_GT_BCF1(ALLELE) ((( (ALLELE) >> 1) - 1) << 1) | (ALLELE & 1)

#define YON_GT_UN_NONE       0       // nothing
#define YON_GT_UN_INT        1       // rcds
#define YON_GT_UN_BCF        2|YON_GT_UN_INT // bcf
#define YON_GT_UN_SIMPLE     4|YON_GT_UN_INT // simple
#define YON_GT_UN_BCF_PPA    8|YON_GT_UN_BCF // bcf unpermuted
#define YON_GT_UN_SIMPLE_PPA 16|YON_GT_UN_SIMPLE // bcf unpermuted
#define YON_GT_UN_ALL        (YON_GT_UN_BCF_PPA|YON_GT_UN_SIMPLE_PPA) // everything

// 0 for missing and 1 for sentinel node. Note that the
// sentinel node never occurs in this encoding type.
const uint8_t YON_GT_RLE_RECODE[3] = {2, 3, 0};

// Vcf:INFO names for fields annotated when triggering
// annotation of genotypes.
const std::vector< std::string > YON_GT_ANNOTATE_FIELDS = {"NM","NPM","AN","HWE_P","AC","AF","AC_P","FS_A","F_PIC","HET","MULTI_ALLELIC"};

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

// Primary generic Tachyon FORMAT:GT structure.
struct yon_gt_rcd {
	yon_gt_rcd() : run_length(0), allele(nullptr){}
	~yon_gt_rcd(){ delete [] this->allele; }
	yon_gt_rcd(const yon_gt_rcd& other) = delete; // disallow copy ctor
	yon_gt_rcd& operator=(const yon_gt_rcd& other) = delete; // disallow move assign
	yon_gt_rcd(yon_gt_rcd&& other) : run_length(other.run_length), allele(other.allele){
		other.allele = nullptr;
	}
	yon_gt_rcd& operator=(yon_gt_rcd&& other){
		if(this == &other) return(*this);
		delete this->allele;
		this->allele = other.allele;
		other.allele = nullptr;
		this->run_length = other.run_length;
		return(*this);
	}

	// Rule of 5

	io::BasicBuffer& PrintVcf(io::BasicBuffer& buffer, const uint8_t& n_ploidy){
		if(this->allele[0] == 1){
			buffer += '.';
			return(buffer);
		}
		if(this->allele[0] == 0) buffer += '.';
		else buffer.AddReadble(((this->allele[0] >> 1) - 2));

		for(U32 i = 1; i < n_ploidy; ++i){
			if(this->allele[i] == 1) break;
			buffer += ((this->allele[i] & 1) ? '|' : '/');
			if(this->allele[i] == 0) buffer += '.';
			else buffer.AddReadble(((this->allele[i] >> 1) - 2));
		}
		return(buffer);
	}

	uint32_t run_length;
	//uint8_t alleles; // Determined from base ploidy
	uint8_t* allele; // contains phase at first bit
};

// Forward declare.
struct yon_gt_summary;

struct yon_gt {
    uint8_t  add : 7,
             global_phase : 1;
    uint8_t  shift;
    uint8_t  p, m, method; // bytes per entry, base ploidy, base method
    uint32_t n_s, n_i, n_o;     // number samples, number of entries
    uint8_t  n_allele;
    yon_gt_ppa* ppa; // pointer to ppa
    std::vector< std::vector<uint32_t> >* occ; // pointer to occ
    uint8_t* data; // pointer to data
    uint8_t* d_bcf; // lazy evaluated as Bcf entries (length = base_ploidy * n_samples * sizeof(uint8_t))
    uint8_t* d_bcf_ppa; // lazy evaluation of unpermuted bcf records
    yon_gt_rcd** d_exp; // lazy evaluated from ppa/normal to internal offset (length = n_samples). This can be very expensive if evaluated internally for every record.
    yon_gt_rcd* rcds; // lazy interpreted internal records
    uint32_t* n_occ;
    yon_gt_rcd** d_occ; // lazy evaluation of occ table
    algorithm::IntervalTree<uint32_t, yon_gt_rcd*>* itree; // interval tree for consecutive ranges
    bool dirty;

    typedef yonRawIterator<yon_gt_rcd>       iterator;
	typedef yonRawIterator<const yon_gt_rcd> const_iterator;

    yon_gt() : add(0), global_phase(0), shift(0), p(0), m(0), method(0), n_s(0), n_i(0), n_o(0),
               n_allele(0), ppa(nullptr), occ(nullptr), data(nullptr), d_bcf(nullptr),
			   d_bcf_ppa(nullptr), d_exp(nullptr), rcds(nullptr), n_occ(nullptr), d_occ(nullptr),
			   itree(nullptr), dirty(false)
    {}

    ~yon_gt();

    bool Evaluate(void){
    	if(this->method == 1) return(this->EvaluateRecordsM1());
    	else if(this->method == 2) return(this->EvaluateRecordsM2());
    	else if(this->method == 4) return(this->EvaluateRecordsM4());
    	else {
    		std::cerr << "not implemented method " << (int)this->method << std::endl;
    	}
    	return false;
    }

    bool EvaluateRecordsM1(){
    	switch(this->p){
    	case(1): return(this->EvaluateRecordsM1_<uint8_t>());
    	case(2): return(this->EvaluateRecordsM1_<uint16_t>());
    	case(4): return(this->EvaluateRecordsM1_<uint32_t>());
    	case(8): return(this->EvaluateRecordsM1_<uint64_t>());
    	default:
    		std::cerr << "illegal primitive in EvaluateRecordsM1" << std::endl;
    		exit(1);
    	}
    }

    template <class T>
    bool EvaluateRecordsM1_(){
    	if(this->rcds != nullptr) delete [] this->rcds;
    	assert(this->m == 2);
    	assert(this->n_allele == 2);

    	// Allocate memory for new records.
    	this->rcds = new yon_gt_rcd[this->n_i];

    	// Reinterpret byte stream into the appropriate actual
    	// primitive type.
    	const T* r_data = reinterpret_cast<const T*>(this->data);

    	// Keep track of the cumulative number of genotypes observed
    	// as a means of asserting correctness.
    	uint64_t n_total = 0;

    	// Iterate over the internal run-length encoded genotypes
    	// and populate the rcds structure.
		for(uint32_t i = 0; i < this->n_i; ++i){
			BYTE phasing = 0;
			if(add) phasing = r_data[i] & 1;
			else    phasing = this->global_phase;

			this->rcds[i].run_length = YON_GT_RLE_LENGTH(r_data[i], shift, add);
			this->rcds[i].allele = new uint8_t[2];
			this->rcds[i].allele[0] = YON_GT_RLE_ALLELE_B(r_data[i], shift, add);
			this->rcds[i].allele[1] = YON_GT_RLE_ALLELE_A(r_data[i], shift, add);
			// Store an allele encoded as (ALLELE << 1 | phasing).
			this->rcds[i].allele[0] = (YON_GT_RLE_RECODE[this->rcds[i].allele[0]] << 1) | phasing;
			this->rcds[i].allele[1] = (YON_GT_RLE_RECODE[this->rcds[i].allele[1]] << 1) | phasing;
			n_total += this->rcds[i].run_length;
		}
		assert(n_total == this->n_s);
    }

    bool EvaluateRecordsM2(){
		switch(this->p){
		case(1): return(this->EvaluateRecordsM2_<uint8_t>());
		case(2): return(this->EvaluateRecordsM2_<uint16_t>());
		case(4): return(this->EvaluateRecordsM2_<uint32_t>());
		case(8): return(this->EvaluateRecordsM2_<uint64_t>());
		default:
			std::cerr << "illegal primitive in EvaluateRecordsM2" << std::endl;
			exit(1);
		}
	}

	template <class T>
	bool EvaluateRecordsM2_(){
		if(this->rcds != nullptr) delete [] this->rcds;
		assert(this->m == 2);

		// Allocate memory for new records.
		this->rcds = new yon_gt_rcd[this->n_i];

		// Reinterpret byte stream into the appropriate actual
		// primitive type.
		const T* r_data = reinterpret_cast<const T*>(this->data);

		// Keep track of the cumulative number of genotypes observed
		// as a means of asserting correctness.
		uint64_t n_total = 0;

		// Iterate over the internal run-length encoded genotypes
		// and populate the rcds structure.
		for(uint32_t i = 0; i < this->n_i; ++i){
			BYTE phasing = 0;
			if(add) phasing = r_data[i] & 1;
			else    phasing = this->global_phase;

			this->rcds[i].run_length = YON_GT_RLE_LENGTH(r_data[i], shift, add);
			this->rcds[i].allele = new uint8_t[2];
			this->rcds[i].allele[0] = YON_GT_RLE_ALLELE_B(r_data[i], shift, add);
			this->rcds[i].allele[1] = YON_GT_RLE_ALLELE_A(r_data[i], shift, add);
			// Store an allele encoded as (ALLELE << 1 | phasing).
			this->rcds[i].allele[0] = (this->rcds[i].allele[0] << 1) | phasing;
			this->rcds[i].allele[1] = (this->rcds[i].allele[1] << 1) | phasing;
			n_total += this->rcds[i].run_length;
		}
		assert(n_total == this->n_s);
	}

	 bool EvaluateRecordsM4(){
		switch(this->p){
		case(1): return(this->EvaluateRecordsM4_<uint8_t>());
		case(2): return(this->EvaluateRecordsM4_<uint16_t>());
		case(4): return(this->EvaluateRecordsM4_<uint32_t>());
		case(8): return(this->EvaluateRecordsM4_<uint64_t>());
		default:
			std::cerr << "illegal primitive in EvaluateRecordsM1" << std::endl;
			exit(1);
		}
	}

	template <class T>
	bool EvaluateRecordsM4_(){
		if(this->rcds != nullptr) delete [] this->rcds;
		assert(this->m != 2);

		// Allocate memory for new records.
		this->rcds = new yon_gt_rcd[this->n_i];

		// Keep track of the cumulative number of genotypes observed
		// as a means of asserting correctness.
		uint64_t n_total  = 0;
		uint64_t b_offset = 0;

		// Iterate over the internal run-length encoded genotypes
		// and populate the rcds structure.
		for(uint32_t i = 0; i < this->n_i; ++i){
			T* run_length = reinterpret_cast<T*>(&this->data[b_offset]);
			b_offset += sizeof(T);

			this->rcds[i].run_length = *run_length;
			this->rcds[i].allele = new uint8_t[this->m];
			for(U32 j = 0; j < this->m; ++j, ++b_offset)
				this->rcds[i].allele[j] = this->data[b_offset];

			n_total += this->rcds[i].run_length;
		}
		assert(n_total == this->n_s);
	}

	// Requires base evaluation to rcds structures first.
	bool EvaluateBcf(){
		if(this->rcds == nullptr){
			std::cerr << "have to evaluate rcds first" << std::endl;
			return false;
		}

		switch(this->p){
		case(1): return(this->EvaluateBcf_<uint8_t>());
		case(2): return(this->EvaluateBcf_<uint16_t>());
		case(4): return(this->EvaluateBcf_<uint32_t>());
		case(8): return(this->EvaluateBcf_<uint64_t>());
		default:
			std::cerr << "illegal primitive in EvaluateBcf" << std::endl;
			exit(1);
		}
	}

	template <class T>
	bool EvaluateBcf_(){
		assert(this->rcds != nullptr);
		if(this->d_bcf != nullptr) delete [] this->d_bcf;

		uint64_t cum_pos = 0;
		this->d_bcf = new uint8_t[this->m * this->n_s];
		for(uint32_t i = 0; i < this->n_i; ++i){
			for(uint32_t j = 0; j < this->rcds[i].run_length; ++j){
				for(uint32_t k = 0; k < this->m; ++k, ++cum_pos){
					this->d_bcf[cum_pos] = this->rcds[i].allele[k]; // Todo: recode back from YON-style to htslib style (e.g. 0,1 special meaning in YON)
				}
			}
		}
		assert(cum_pos == this->n_s * this->m);
		return true;
	}

	bool EvaluatePpaFromBcf(const uint8_t* d_bcf_pre){
		assert(d_bcf_pre != nullptr);
		assert(this->ppa != nullptr);

		if(this->d_bcf_ppa != nullptr) delete [] this->d_bcf_ppa;
		this->d_bcf_ppa = new uint8_t[this->n_s * this->m];

		uint64_t cum_pos = 0;
		for(uint32_t i = 0; i < this->n_s; ++i){
			const uint32_t& ppa_target = this->ppa->ordering[i];
			const uint8_t* bcf_target = &d_bcf_pre[ppa_target * this->m];

			for(uint32_t k = 0; k < this->m; ++k, ++cum_pos){
				this->d_bcf_ppa[cum_pos] = bcf_target[k];
			}
		}

		assert(cum_pos == this->n_s * this->m);
		return true;
	}

	bool EvaluateIntervalTree(void){
		assert(this->rcds != nullptr);
		if(this->itree != nullptr) delete this->itree;

		std::vector< algorithm::Interval<uint32_t, yon_gt_rcd*> > intervals;
		uint64_t cum_pos = 0;
		for(uint32_t i = 0; i < this->n_i; ++i){
			intervals.push_back(algorithm::Interval<uint32_t, yon_gt_rcd*>(
					cum_pos,
					cum_pos + this->rcds[i].run_length,
					&this->rcds[i])
				);
			cum_pos += this->rcds[i].run_length;
		}
		this->itree = new algorithm::IntervalTree<uint32_t, yon_gt_rcd*>(std::move(intervals));

		assert(cum_pos == this->n_s);
		return true;
	}

	/**<
	 * Lazy evaluated (expands) (possibly) run-length encoded rcds structures
	 * into a vector of pointers corresponding to distinct samples in the order
	 * they were stored. This unpermutation (restoration) of order requires the
	 * ppa array (stored at YON_BLK_PPA).
	 *
	 * This function should be considered private as the generic wrapper function
	 * Expand() calls this function if ppa is available or ExpandRecords()
	 * otherwise.
	 * @return Always returns TRUE if the critical assertions passes.
	 */
	bool ExpandRecordsPpa(void){
		assert(this->rcds != nullptr);
		assert(this->ppa != nullptr);

		if(this->d_exp != nullptr) delete [] this->d_exp;
		this->d_exp = new yon_gt_rcd*[this->n_s];

		uint64_t cum_sample = 0;
		uint64_t cum_offset = 0;
		for(uint32_t i = 0; i < this->n_i; ++i){
			for(uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample){
				const uint32_t& target_ppa = this->ppa->at(cum_sample);
				this->d_exp[target_ppa] = &this->rcds[i];
				cum_offset += this->p * this->m;
			}
		}
		assert(cum_sample == this->n_s);
		assert(cum_offset == this->n_s * this->m * this->p);
		return true;
	}

	bool ExpandRecordsPpaExternal(yon_gt_rcd** d_expe){
		assert(this->rcds != nullptr);
		assert(this->ppa != nullptr);

		uint64_t cum_sample = 0;
		uint64_t cum_offset = 0;
		for(uint32_t i = 0; i < this->n_i; ++i){
			for(uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample){
				const uint32_t& target_ppa = this->ppa->at(cum_sample);
				d_expe[target_ppa] = &this->rcds[i];
				cum_offset += this->p * this->m;
			}
		}
		assert(cum_sample == this->n_s);
		assert(cum_offset == this->n_s * this->m * this->p);
		return true;
	}

	bool Expand(void){
		if(this->ppa != nullptr)
			return(this->ExpandRecordsPpa());
		else return(this->ExpandRecords());
	}

	bool ExpandExternal(yon_gt_rcd** d_expe){
		if(this->ppa != nullptr)
			return(this->ExpandRecordsPpaExternal(d_expe));
		else return(this->ExpandRecordsExternal(d_expe));
	}


	bool ExpandRecords(void){
		assert(this->rcds != nullptr);

		if(this->d_exp != nullptr) delete [] this->d_exp;
		this->d_exp = new yon_gt_rcd*[this->n_s];

		uint64_t cum_sample = 0;
		uint64_t cum_offset = 0;
		for(uint32_t i = 0; i < this->n_i; ++i){
			for(uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample){
				this->d_exp[cum_sample] = &this->rcds[i];
				cum_offset += this->p * this->m;
			}
		}
		assert(cum_sample == this->n_s);
		assert(cum_offset == this->n_s * this->m * this->p);
		return true;
	}

	bool ExpandRecordsExternal(yon_gt_rcd** d_expe){
		assert(this->rcds != nullptr);

		uint64_t cum_sample = 0;
		uint64_t cum_offset = 0;
		for(uint32_t i = 0; i < this->n_i; ++i){
			for(uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample){
				d_expe[cum_sample] = &this->rcds[i];
				cum_offset += this->p * this->m;
			}
		}
		assert(cum_sample == this->n_s);
		assert(cum_offset == this->n_s * this->m * this->p);
		return true;
	}

    template <class T>
    inline T& GetPrimitive(const uint32_t sample){ return(*reinterpret_cast<T*>(&this->data[sample])); }

    template <class T>
	inline T& GetPrimitivePpa(const uint32_t sample){ return(this->ppa[*reinterpret_cast<T*>(&this->data[sample])]); }

    // Iterator
   inline iterator begin(){ return iterator(&this->rcds[0]); }
   inline iterator end()  { return iterator(&this->rcds[this->n_i]); }
   inline const_iterator begin()  const{ return const_iterator(&this->rcds[0]); }
   inline const_iterator end()    const{ return const_iterator(&this->rcds[this->n_i]); }
   inline const_iterator cbegin() const{ return const_iterator(&this->rcds[0]); }
   inline const_iterator cend()   const{ return const_iterator(&this->rcds[this->n_i]); }

   /**<
    * Transform lazy-evaluated Tachyon genotype encodings (d_exp)
    * to htslib bcf1_t genotype encodings.
    * @param rec Input bcf1_t record.
    * @param hdr Input htslib bcf header.
    * @return    Returns the pointer to the input bcf1_t record.
    */
   bcf1_t* UpdateHtslibGenotypes(bcf1_t* rec, bcf_hdr_t* hdr) const{
	   assert(this->d_exp != nullptr);

	   int32_t* tmpi = new int32_t[this->n_s*this->m];
	   uint32_t gt_offset = 0;
	   for(U32 i = 0; i < this->n_s; ++i){
		   for(U32 j = 0; j < this->m; ++j, ++gt_offset){
			   if(this->d_exp[i]->allele[j] == 0)      tmpi[gt_offset] = 0;
			   else if(this->d_exp[i]->allele[j] == 1) tmpi[gt_offset] = 1;
			   else tmpi[gt_offset] = YON_GT_BCF1(this->d_exp[i]->allele[j]);
		   }
	   }
	   assert(gt_offset == this->n_s*this->m);

	   bcf_update_genotypes(hdr, rec, tmpi, this->n_s*this->m);
	   delete [] tmpi;
	   return(rec);
   }
};

struct yon_gt_summary_obj{
	yon_gt_summary_obj() : n_cnt(0), children(nullptr){}
	~yon_gt_summary_obj(){ delete [] this->children; }

	inline yon_gt_summary_obj& operator[](const uint32_t pos){ return(this->children[pos]); }
	inline const yon_gt_summary_obj& operator[](const uint32_t pos) const{ return(this->children[pos]); }

	uint64_t n_cnt;
	yon_gt_summary_obj* children;
};

struct yon_gt_summary_rcd {
	yon_gt_summary_rcd() :
		n_ploidy(0), n_ac_af(0), n_fs(0), ac(nullptr), af(nullptr),
		nm(0), npm(0), an(0), ac_p(nullptr), fs_a(nullptr), hwe_p(0),
		f_pic(0), heterozygosity(0)
	{}
	~yon_gt_summary_rcd(){
		delete [] this->ac; delete [] this->af;
		delete [] this->fs_a;
		// Do not delete ac_p as it is borrowed
	}

	bcf1_t* UpdateHtslibVcfRecord(bcf1_t* rec, bcf_hdr_t* hdr) const{
		utility::UpdateHtslibVcfRecordInfo(rec, hdr, "NM", (const U64*)&this->nm, 1);
		utility::UpdateHtslibVcfRecordInfo(rec, hdr, "NPM", (const U64*)&this->npm, 1);
		utility::UpdateHtslibVcfRecordInfo(rec, hdr, "AN", (const U64*)&this->an, 1);
		utility::UpdateHtslibVcfRecordInfo(rec, hdr, "HWE_P", (const double*)&this->hwe_p, 1);

		if(this->n_ac_af > 4)
			bcf_update_info_flag(hdr, rec, "MULTI_ALLELIC", NULL, 1);

		if(this->n_ac_af > 2){
			utility::UpdateHtslibVcfRecordInfo(rec, hdr, "AC", (const U64*)&this->ac[2], this->n_ac_af - 2);
			utility::UpdateHtslibVcfRecordInfo(rec, hdr, "AF", (const double*)&this->af[2], this->n_ac_af - 2);
		}

		std::vector<uint64_t> ac_p;
		for(U32 p = 0; p < this->n_ploidy; ++p){
			for(U32 i = 2; i < this->n_ac_af; ++i){
				ac_p.push_back(this->ac_p[p][i]);
			}
		}
		utility::UpdateHtslibVcfRecordInfo(rec, hdr, "AC_P", (const U64*)ac_p.data(), ac_p.size());

		if(this->fs_a != nullptr)
			utility::UpdateHtslibVcfRecordInfo(rec, hdr, "FS_A", this->fs_a, this->n_fs);

		if(this->n_ploidy == 2){
			utility::UpdateHtslibVcfRecordInfo(rec, hdr, "FPIC", &this->f_pic, 1);
			utility::UpdateHtslibVcfRecordInfo(rec, hdr, "HET", &this->heterozygosity, 1);
		}

		return(rec);
	}

	io::BasicBuffer& PrintVcf(io::BasicBuffer& buffer){
		buffer +=  "NM=";     buffer.AddReadble((U64)this->nm);
		buffer += ";NPM=";    buffer.AddReadble((U64)this->npm);
		buffer += ";AN=";     buffer.AddReadble((U64)this->an);
		buffer += ";HWE_P=";  buffer.AddReadble((double)this->hwe_p);
		if(this->n_ac_af > 4) buffer += ";MULTI_ALLELIC";

		if(this->n_ac_af > 2){
			buffer += ";AC="; buffer.AddReadble((U64)this->ac[2]);
			for(U32 i = 3; i < this->n_ac_af; ++i){
				buffer += ','; buffer.AddReadble((U64)this->ac[i]);
			}
			buffer += ";AF="; buffer.AddReadble((double)this->af[2]);
			for(U32 i = 3; i < this->n_ac_af; ++i){
				buffer += ','; buffer.AddReadble((double)this->af[i]);
			}
		}

		buffer += ";AC_P=";
		buffer.AddReadble((U64)this->ac_p[0][2]);
		for(U32 i = 3; i < this->n_ac_af; ++i){
			buffer += ','; buffer.AddReadble((U64)this->ac_p[0][i]);
		}

		for(U32 p = 1; p < this->n_ploidy; ++p){
			for(U32 i = 2; i < this->n_ac_af; ++i){
				buffer += ','; buffer.AddReadble((U64)this->ac_p[p][i]);
			}
		}

		if(this->fs_a != nullptr){
			buffer += ";FS_A=";
			buffer.AddReadble((double)this->fs_a[0]);
			for(U32 i = 1; i < this->n_fs; ++i){
				buffer += ','; buffer.AddReadble((double)this->fs_a[i]);
			}
		}

		if(this->n_ploidy == 2){
			buffer += ";F_PIC=";
			buffer.AddReadble((double)this->f_pic);
			buffer += ";HET=";
			buffer.AddReadble((double)this->heterozygosity);
		}

		return(buffer);
	}

	uint8_t n_ploidy;
	uint32_t n_ac_af;
	uint32_t n_fs;
	uint64_t* ac; // allele counts
	double* af; // allele frequency
	uint64_t nm; // number of non-sentinel, non-missing symbols
	uint64_t npm;// number of missing symbols
	uint64_t an; // number of non-sentinel symbols
	// ac_p is borrowed from summary
	uint64_t** ac_p;
	double* fs_a;// fisher strand test p
	double hwe_p;// hardy-weinberg p
	double f_pic;
	double heterozygosity;
};

struct yon_gt_summary{
	yon_gt_summary(void) :
		n_ploidy(0),
		n_alleles(0),
		alleles(nullptr),
		alleles_strand(nullptr),
		gt(nullptr),
		d(nullptr)
	{

	}

	yon_gt_summary(const uint8_t base_ploidy, const uint8_t n_alleles) :
		n_ploidy(base_ploidy),
		n_alleles(n_alleles + 2),
		alleles(new uint64_t[this->n_alleles]),
		alleles_strand(new uint64_t*[this->n_ploidy]),
		gt(new yon_gt_summary_obj[this->n_alleles]),
		d(nullptr)
	{
		memset(this->alleles, 0, sizeof(uint64_t)*this->n_alleles);
		for(U32 i = 0; i < this->n_ploidy; ++i){
			this->alleles_strand[i] = new uint64_t[this->n_alleles];
			memset(this->alleles_strand[i], 0, sizeof(uint64_t)*this->n_alleles);
		}

		// Add layers to the root node.
		for(U32 i = 0; i < this->n_alleles; ++i)
			this->AddGenotypeLayer(&this->gt[i], 1);
	}

	~yon_gt_summary(){
		delete [] this->alleles;
		delete [] this->gt;
		if(this->alleles_strand != nullptr){
			for(U32 i = 0; i < this->n_ploidy; ++i)
				delete this->alleles_strand[i];
			delete [] this->alleles_strand;
		}
		delete this->d;
	}

	/**<
	 * Recursively add layers to the full trie in order to represent
	 * a possible ploidy-dimensional matrix at the leafs. The memory
	 * cost of the trie is O(n_alleles ^ ploidy) and allows the prefix
	 * lookup of any genotype.
	 *
	 * @param target
	 * @param depth
	 * @return
	 */
	bool AddGenotypeLayer(yon_gt_summary_obj* target, uint8_t depth);

	/**<
	 *
	 * @param gt
	 * @return
	 */
	yon_gt_summary& operator+=(const yon_gt& gt);

	// Accessors to internal data.
	inline uint64_t* GetAlleleCountsRaw(void){ return(this->alleles); }
	inline const uint64_t* GetAlleleCountsRaw(void) const{ return(this->alleles); }

	std::vector<uint64_t> GetAlleleCounts(void) const;
	std::vector< std::pair<uint64_t,double> > GetAlleleCountFrequency(void) const;
	std::vector< std::vector<uint64_t> > GetAlleleStrandCounts(void) const;
	bool GetGenotype(std::vector<uint64_t>& data,
	                 yon_gt_summary_obj* target,
	                 uint8_t depth) const;

	// Todo: unfinished.
	std::vector<yon_gt_rcd> GetGenotypeCounts(bool drop_empty = true) const;
	std::vector<double> GetStrandBiasAlleles(const bool phred_scale = true) const;
	double CalculateHardyWeinberg(void) const;

	/**<
	 *
	 * @return
	 */
	bool LazyEvaluate(void);

	uint8_t    n_ploidy; // base ploidy at site
	uint8_t    n_alleles; // number of alleles
	uint64_t*  alleles; // allele counts
	uint64_t** alleles_strand; // allelic counts per chromosome
	yon_gt_summary_obj* gt; // full genotypic trie with branch-size n_alleles
	yon_gt_summary_rcd* d; // lazy evaluated record
};


}

#endif /* CORE_GENOTYPES_H_ */
