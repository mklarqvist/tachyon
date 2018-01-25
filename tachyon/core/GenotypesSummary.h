#ifndef ENCODERGENOTYPESRLEHELPER_H_
#define ENCODERGENOTYPESRLEHELPER_H_

#include "GTObject.h"

namespace tachyon{
namespace core{

#define DIPLOID_ALLELE_LOOKUP(A,B,shift,mask) ((A & mask) << shift) | (B & mask)
// Forward declare
class GenotypeContainerInterface;
template <class T> class GenotypeContainerDiploidRLE;
template <class T> class GenotypeContainerDiploidSimple;

/**<
 * We cannot decouple the GenotypeContainer objects from
 * this structure as their interpretation is contingent
 * on the internal memory layout of this structure.
 *
 * Note the non-polynomial expansion of `__allele_cells`
 * cost. Number of entries for `n_alleles` {1..10} is shown:
 * 4      16      64     256    1024    4096   16384   65536  262144
 *
 * At `n_alleles` = 14 the memory cost for this object exceeds
 * 2 GB!
 */
struct GenotypeSum{
	typedef GenotypeSum self_type;

	GenotypeSum() :
		n_alleles(5),
		n_cells(pow(2,(this->n_alleles*2))),           // 2^(n_alelles*2) where 2 is ploidy and base is bits
		__genotype_cells(new U64[this->n_cells]),
		__alleleA_cells(new U64[this->n_alleles]),
		__alleleB_cells(new U64[this->n_alleles]),
		__bit_mask((1 << this->n_alleles) - 1),
		__bit_shift_width(this->n_alleles)
	{
		memset(this->__genotype_cells,  0, sizeof(U64)*this->n_cells);
		memset(this->__alleleA_cells,   0, sizeof(U64)*this->n_alleles);
		memset(this->__alleleB_cells,   0, sizeof(U64)*this->n_alleles);
	}

	GenotypeSum(const BYTE& number_of_alleles) :
		n_alleles(number_of_alleles),
		n_cells(pow(2,(this->n_alleles*2))),           // 2^(n_alelles*2) where 2 is ploidy and base is bits
		__genotype_cells(new U64[this->n_cells]),
		__alleleA_cells(new U64[this->n_alleles]),
		__alleleB_cells(new U64[this->n_alleles]),
		__bit_mask((1 << this->n_alleles) - 1),
		__bit_shift_width(this->n_alleles)
	{
		memset(this->__genotype_cells,  0, sizeof(U64)*this->n_cells);
		memset(this->__alleleA_cells,   0, sizeof(U64)*this->n_alleles);
		memset(this->__alleleB_cells,   0, sizeof(U64)*this->n_alleles);
	}

	~GenotypeSum(){
		delete [] this->__genotype_cells;
		delete [] this->__alleleA_cells;
		delete [] this->__alleleB_cells;
	}

	// Overloaded operators
	inline self_type& operator+=(const self_type& other){ // TODO
		// if their sizes are different then we need to map
		// from object 1 to object 2
		// map from SMALLER -> LARGER by expanding bit-width:
		// a) if this is smaller then expand this first to size of other and map 1:1
		// b) if the other object is smaller then map from other to this by bit-expansion
		return(*this);
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidRLE<T>& gt_rle_container){
		const BYTE shift = gt_rle_container.getMeta().isAnyGTMissing()    ? 2 : 1;
		const BYTE add   = gt_rle_container.getMeta().isGTMixedPhasing()  ? 1 : 0;

		// Have to remap to shared type
		// 0 is missing for other types
		// but higher value in RLE
		// remap RLE -> Regular

		for(U32 i = 0; i < gt_rle_container.size(); ++i){
			const BYTE alleleA = (gt_rle_container[i] & ((1 << shift) - 1) << add) >> add;
			const BYTE alleleB = (gt_rle_container[i] & ((1 << shift) - 1) << (add+shift)) >> (add+shift);
			const U32  length  = gt_rle_container[i] >> (2*shift + add);
			this->getGenotypeCell(alleleA, alleleB) += length;
			this->getAlleleA(alleleA) += length;
			this->getAlleleB(alleleB) += length;
		}
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidSimple<T>& gt_simple_container){
		const BYTE shift = ceil(log2(gt_simple_container.getMeta().getNumberAlleles() + 1 + gt_simple_container.getMeta().isAnyGTMissing())); // Bits occupied per allele, 1 value for missing
		const BYTE add   = gt_simple_container.getMeta().isGTMixedPhasing() ? 1 : 0;

		for(U32 i = 0; i < gt_simple_container.size(); ++i){
			const BYTE alleleA = (gt_simple_container[i] & ((1 << shift) - 1) << add) >> add;
			const BYTE alleleB = (gt_simple_container[i] & ((1 << shift) - 1) << (add+shift)) >> (add+shift);
			const U32  length  = gt_simple_container[i] >> (2*shift + add);
			this->getGenotypeCell(alleleA, alleleB) += length;
			this->getAlleleA(alleleA) += length;
			this->getAlleleB(alleleB) += length;
		}
	}

	// Utility
	void resize(const BYTE& n_alleles);
	inline const size_t& size(void) const{ return(this->n_cells); }
	inline void clear(void){
		memset(this->__genotype_cells,  0, sizeof(U64)*this->n_cells);
		memset(this->__alleleA_cells,   0, sizeof(U64)*this->n_alleles);
		memset(this->__alleleB_cells,   0, sizeof(U64)*this->n_alleles);
	}

	// Genotype accessors
	inline U64& operator()(const BYTE& allele1, const BYTE& allele2){
		assert((DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)) < this->n_cells);
		return(this->__genotype_cells[DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)]);
	}

	inline const U64& operator()(const BYTE& allele1, const BYTE& allele2) const{
		assert((DIPLOID_ALLELE_LOOKUP(allele1, allele2, this->__bit_shift_width, this->__bit_mask)) < this->n_cells);
		return(this->__genotype_cells[DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)]);
	}

	inline U64& getGenotypeCell(const BYTE& allele1, const BYTE& allele2){
		assert((DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)) < this->n_cells);
		return(this->__genotype_cells[DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)]);
	}

	inline const U64& getGenotypeCell(const BYTE& allele1, const BYTE& allele2) const{
		assert((DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)) < this->n_cells);
		return(this->__genotype_cells[DIPLOID_ALLELE_LOOKUP(allele1 - 1, allele2 - 1, this->__bit_shift_width, this->__bit_mask)]);
	}

	// Allele accessors
	inline U64& getAlleleA(const BYTE& allele){ return(this->__alleleA_cells[allele]); }
	inline const U64& getAlleleA(const BYTE& allele) const{ return(this->__alleleA_cells[allele]); }
	inline U64& getAlleleB(const BYTE& allele){ return(this->__alleleB_cells[allele]); }
	inline const U64& getAlleleB(const BYTE& allele) const{ return(this->__alleleB_cells[allele]); }
	inline U64 getAllele(const BYTE& allele) const{ return(this->__alleleA_cells[allele] + this->__alleleB_cells[allele]); }

	//
	U64 alleleCount(void) const{
		U64 total = 0;
		for(U32 i = 0; i < this->n_alleles; ++i)
			total += this->getAllele(i);
		return(total);
	}

	U64 genotypeCount(void) const{
		U64 total = 0;
		for(U32 i = 0; i < this->n_cells; ++i)
			total += this->__genotype_cells[i];
		return(total);
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << "alleles = [" << entry.getAllele(0);
		for(U32 i = 1; i < entry.n_alleles; ++i){
			out << ',' << entry.getAllele(i);
		}
		out << "] allelesA = [" << entry.getAllele(0);
		for(U32 i = 1; i < entry.n_alleles; ++i){
			out << ',' << entry.getAlleleA(i);
		}
		out << "] allelesB = [" << entry.getAllele(0);
		for(U32 i = 1; i < entry.n_alleles; ++i){
			out << ',' << entry.getAlleleB(i);
		}

		out << "] genotypes = [";
		for(U32 i = 0; i < entry.n_alleles; ++i){
			for(U32 j = 0; j < entry.n_alleles; ++j){
				const U64& target_cell = entry.getGenotypeCell(i,j);
				if(target_cell)
					out << i << "," << j << ":" << target_cell << ',';
			}
		}
		out << "]";

		return(out);
	}

private:
	BYTE   n_alleles;         // number of alleles (including NA and missing)
	size_t n_cells;           // length of `__genotype_cells`
	U64*   __genotype_cells;  // length is `n_cells`
	U64*   __alleleA_cells;   // length is `n_alleles`
	U64*   __alleleB_cells;   // length is `n_alleles`
	BYTE   __bit_mask;        // guard against undesired residual overflow in packing
	BYTE   __bit_shift_width; // bit shift size for packing
};

// value1 << ceiling(log2(n_alleles))

struct GenotypesSummary{
	typedef GenotypesSummary self_type;

	GenotypesSummary(void) :
		MAF(0),
		MGF(0),
		HWE_P(0),
		missingValues(0),
		phased(false),
		expectedSamples(0)
	{
		memset(&this->countsGenotypes[0], 0, sizeof(U64)*16);
		memset(&this->countsAlleles[0],   0, sizeof(U64)*3);
	}

	GenotypesSummary(const U64 expectedSamples) :
		MAF(0),
		MGF(0),
		HWE_P(0),
		missingValues(0),
		phased(false),
		expectedSamples(expectedSamples)
	{
		memset(&this->countsGenotypes[0], 0, sizeof(U64)*16);
		memset(&this->countsAlleles[0],   0, sizeof(U64)*3);
	}
	~GenotypesSummary(){}

	void setExpected(const U32 expected){ this->expectedSamples = expected; }
	U64& operator[](const U32& p){ return(this->countsGenotypes[p]); }

	inline void reset(){
		this->countsGenotypes[0] = 0;
		this->countsGenotypes[1] = 0;
		this->countsGenotypes[4] = 0;
		this->countsGenotypes[5] = 0;
		this->countsAlleles[0]   = 0; // p
		this->countsAlleles[1]   = 0; // q
		this->countsAlleles[2]   = 0; // missing
	}
	inline const bool& hasMissing(void) const{ return(this->missingValues); }

	double calculateMAF(void){
		if(this->countsAlleles[0] > this->countsAlleles[1])
			return(this->countsAlleles[1]/((double)this->countsAlleles[0]+this->countsAlleles[1]));
		else
			return(this->countsAlleles[0]/((double)this->countsAlleles[0]+this->countsAlleles[1]));
	}

	void calculateMGF(void){
		// Find the largest non-missing value
		U64 curMax = 0;
		U64* target; // To compare pointer address
		if(this->countsGenotypes[0] > curMax){
			curMax = this->countsGenotypes[0];
			target = &this->countsGenotypes[0];
		}
		if(this->countsGenotypes[1] > curMax){
			curMax = this->countsGenotypes[1];
			target = &this->countsGenotypes[1];
		}
		if(this->countsGenotypes[4] > curMax){
			curMax = this->countsGenotypes[4];
			target = &this->countsGenotypes[4];
		}
		if(this->countsGenotypes[5] > curMax){
			curMax = this->countsGenotypes[5];
			target = &this->countsGenotypes[5];
		}

		// Find next largest non-missing value
		U64 curMax2 = 0;
		if(this->countsGenotypes[0] >= curMax2 && &this->countsGenotypes[0] != target)
			curMax2 = this->countsGenotypes[0];
		if(this->countsGenotypes[1] >= curMax2 && &this->countsGenotypes[1] != target)
			curMax2 = this->countsGenotypes[1];
		if(this->countsGenotypes[4] >= curMax2 && &this->countsGenotypes[4] != target)
			curMax2 = this->countsGenotypes[4];
		if(this->countsGenotypes[5] >= curMax2 && &this->countsGenotypes[5] != target)
			curMax2 = this->countsGenotypes[5];

		this->MGF = (double)curMax2/this->expectedSamples;
	}

	bool determinePhase(const char& separator){
		if(separator == '/'){
			this->phased = false;
			return true;
		} else if(separator == '|'){
			this->phased = true;
			return true;
		} else return false;
	}

	/*
	// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
	//
	// Written by Jan Wigginton
	// Modified to use Tachyon data
	*/
	void calculateHardyWeinberg(void){
		U64 obs_hets = this->countsGenotypes[1] + this->countsGenotypes[4];
		U64 obs_hom1 = this->countsGenotypes[0];
		U64 obs_hom2 = this->countsGenotypes[5];

		U64 obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
		U64 obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

		int64_t rare_copies = 2 * obs_homr + obs_hets;
		int64_t genotypes   = obs_hets + obs_homc + obs_homr;

		double* het_probs = new double[rare_copies + 1];

		int64_t i;
		for (i = 0; i <= rare_copies; ++i)
			het_probs[i] = 0.0;

		/* start at midpoint */
		int64_t mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

		/* check to ensure that midpoint and rare alleles have same parity */
		if ((rare_copies & 1) ^ (mid & 1))
			++mid;

		int64_t curr_hets = mid;
		int64_t curr_homr = (rare_copies - mid) / 2;
		int64_t curr_homc = genotypes - curr_hets - curr_homr;

		het_probs[mid] = 1.0;
		double sum = het_probs[mid];
		for (curr_hets = mid; curr_hets > 1; curr_hets -= 2){
			het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
							   / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
			sum += het_probs[curr_hets - 2];

			/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
			++curr_homr;
			++curr_homc;
		}

		curr_hets = mid;
		curr_homr = (rare_copies - mid) / 2;
		curr_homc = genotypes - curr_hets - curr_homr;
		for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2){
			het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
							/((curr_hets + 2.0) * (curr_hets + 1.0));
			sum += het_probs[curr_hets + 2];

			/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
			--curr_homr;
			--curr_homc;
		}

		for (i = 0; i <= rare_copies; i++)
			het_probs[i] /= sum;

		double p_hwe = 0.0;
		/*  p-value calculation for p_hwe  */
		for (i = 0; i <= rare_copies; i++){
			if (het_probs[i] > het_probs[obs_hets])
				continue;

			p_hwe += het_probs[i];
		}

		p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

		delete [] het_probs;

		this->HWE_P = p_hwe;
	}

	friend std::ostream& operator<<(std::ostream& os, const self_type& self){
		os << self.countsGenotypes[0] << '\t' << self.countsGenotypes[1] << '\t' << self.countsGenotypes[4] << '\t' << self.countsGenotypes[5] << '\t' << self.MGF << '\t' << self.HWE_P;
		return(os);
	}

	inline const U64 countAlleles(void) const{ return(this->countsAlleles[0] + this->countsAlleles[1] + this->countsAlleles[2]); }

	U64 countsGenotypes[16];
	U64 countsAlleles[3];
	float MAF;
	float MGF;
	float HWE_P;
	bool missingValues;
	bool phased;
	U64 expectedSamples;
};

}
}



#endif /* ENCODERGENOTYPESRLEHELPER_H_ */
