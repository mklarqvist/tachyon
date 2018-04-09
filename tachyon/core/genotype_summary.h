#ifndef ENCODERGENOTYPESRLEHELPER_H_
#define ENCODERGENOTYPESRLEHELPER_H_

#include "genotype_object.h"

namespace tachyon{
namespace containers{

// Forward declare
class GenotypeContainerInterface;
template <class T> class GenotypeContainerDiploidRLE;
template <class T> class GenotypeContainerDiploidSimple;
template <class T> class GenotypeContainerDiploidBCF;

// Remaps
const BYTE TACHYON_GT_SUMMARY_REMAP[4] = {2, 3, 1, 1}; // 0 = EOV does not exist in this encoding

struct GTSummaryObject{
	GTSummaryObject() : counts_(0), countsA_(0), countsB_(0){}
	~GTSummaryObject(){}

	void reset(void){
		this->counts_ = 0;
		this->countsA_ = 0;
		this->countsB_ = 0;
	}

	void operator+=(const U64& value){ this->counts_ += value; }

	U64 counts_;
	U64 countsA_;
	U64 countsB_;
};

struct GTSummary{
	typedef U64 value_type;

public:
	GTSummary() :
		n_alleles_(7),
		matrix_(new value_type*[this->n_alleles_]),
		vectorA_(new value_type[this->n_alleles_]),
		vectorB_(new value_type[this->n_alleles_])
	{
		for(U32 i = 0; i < this->n_alleles_; ++i){
			this->matrix_[i] = new value_type[this->n_alleles_];
			memset(this->matrix_[i], 0, sizeof(value_type)*this->n_alleles_);
		}
		memset(this->vectorA_, 0, sizeof(value_type)*this->n_alleles_);
		memset(this->vectorB_, 0, sizeof(value_type)*this->n_alleles_);
	}

	GTSummary(const BYTE n_alleles) :
		n_alleles_(n_alleles + 2),
		matrix_(new value_type*[this->n_alleles_]),
		vectorA_(new value_type[this->n_alleles_]),
		vectorB_(new value_type[this->n_alleles_])
	{
		for(U32 i = 0; i < this->n_alleles_; ++i){
			this->matrix_[i] = new value_type[this->n_alleles_];
			memset(this->matrix_[i], 0, sizeof(value_type)*this->n_alleles_);
		}
		memset(this->vectorA_, 0, sizeof(value_type)*this->n_alleles_);
		memset(this->vectorB_, 0, sizeof(value_type)*this->n_alleles_);
	}

	~GTSummary(){
		for(U32 i = 0; i < this->n_alleles_; ++i)
			delete [] this->matrix_[i];
		delete [] this->matrix_;
		delete [] this->vectorA_;
		delete [] this->vectorB_;
	}

	void clear(void){
		for(U32 i = 0; i < this->n_alleles_; ++i){
			memset(this->matrix_[i], 0, sizeof(value_type)*this->n_alleles_);
		}
		memset(this->vectorA_, 0, sizeof(value_type)*this->n_alleles_);
		memset(this->vectorB_, 0, sizeof(value_type)*this->n_alleles_);
	}

	void printDiploid(std::ostream& stream){
		stream << this->matrix_[0][0];
		for(U32 j = 1; j < this->n_alleles_; ++j){
			if(this->matrix_[0][j]) stream << '\t' << this->matrix_[0][j];
		}

		for(U32 i = 1; i < this->n_alleles_; ++i){
			for(U32 j = 0; j < this->n_alleles_; ++j){
				if(this->matrix_[i][j]) stream << '\t' << this->matrix_[i][j];
			}
		}
	}

	U64 alleleCount(void) const{
		U64 total = 0;
		for(U32 i = 2; i < this->n_alleles_; ++i){
			total += this->vectorA_[i];
			total += this->vectorB_[i];
		}
		return(total);
	}

	U64 genotypeCount(void) const{
		U64 total = 0;
		for(U32 i = 2; i < this->n_alleles_; ++i){
			for(U32 j = 2; j < this->n_alleles_; ++j){
				total += this->matrix_[i][j];
			}
		}
		return(total);
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidRLE<T>& gt_rle_container){
		const BYTE shift = gt_rle_container.getMeta().isAnyGTMissing()    ? 2 : 1;
		const BYTE add   = gt_rle_container.getMeta().isGTMixedPhasing()  ? 1 : 0;

		for(U32 i = 0; i < gt_rle_container.size(); ++i){
			const U64 length  = YON_GT_RLE_LENGTH(gt_rle_container.at(i), shift, add);
			const BYTE alleleA = YON_GT_RLE_ALLELE_A(gt_rle_container.at(i), shift, add);
			const BYTE alleleB = YON_GT_RLE_ALLELE_B(gt_rle_container.at(i), shift, add);
			this->matrix_[TACHYON_GT_SUMMARY_REMAP[alleleA]][TACHYON_GT_SUMMARY_REMAP[alleleB]] += length;
			this->vectorA_[TACHYON_GT_SUMMARY_REMAP[alleleA]] += length;
			this->vectorB_[TACHYON_GT_SUMMARY_REMAP[alleleB]] += length;
		}
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidSimple<T>& gt_simple_container){
		const BYTE shift = ceil(log2(gt_simple_container.getMeta().getNumberAlleles() + 1 + gt_simple_container.getMeta().isAnyGTMissing() + gt_simple_container.getMeta().isMixedPloidy())); // Bits occupied per allele, 1 value for missing
		const BYTE add   = gt_simple_container.getMeta().isGTMixedPhasing() ? 1 : 0;
		const BYTE matrix_add = !gt_simple_container.getMeta().isMixedPloidy();
		if(gt_simple_container.getMeta().n_alleles + 2 > this->n_alleles_){
			std::cerr << "too many alleles: " << gt_simple_container.getMeta().n_alleles + 2 << "/" << (int)this->n_alleles_ << std::endl;
			return;
		}

		for(U32 i = 0; i < gt_simple_container.size(); ++i){
			const U64 length = YON_GT_RLE_LENGTH(gt_simple_container.at(i), shift, add);
			BYTE alleleA    = YON_GT_RLE_ALLELE_A(gt_simple_container.at(i), shift, add);
			BYTE alleleB    = YON_GT_RLE_ALLELE_B(gt_simple_container.at(i), shift, add);

			//std::cerr << "adding: " << (int)alleleA << "+" << (alleleA > 0 ? matrix_add : 0) << "/" << (int)alleleB << "+" << (alleleB > 0 ? matrix_add : 0) << ": " << (int)matrix_add << std::endl;

			alleleA += alleleA > 0 ? matrix_add : 0;
			alleleB += alleleB > 0 ? matrix_add : 0;

			this->matrix_[alleleA][alleleB] += length;
			this->vectorA_[alleleA] += length;
			this->vectorB_[alleleB] += length;
		}
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidBCF<T>& gt_diploid_bcf_container){

	}

public:
	BYTE  n_alleles_; // number of alleles (including EOV and missing)
	value_type** matrix_;
	value_type*  vectorA_;
	value_type*  vectorB_;
};

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

	/**<
	// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
	//
	// Written by Jan Wigginton
	// Modified to use Tachyon data by Marcus D. R. Klarqvist
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
		for (i = 0; i <= rare_copies; ++i) het_probs[i] = 0.0;

		// start at midpoint
		int64_t mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

		// check to ensure that midpoint and rare alleles have same parity
		if ((rare_copies & 1) ^ (mid & 1)) ++mid;

		int64_t curr_hets = mid;
		int64_t curr_homr = (rare_copies - mid) / 2;
		int64_t curr_homc = genotypes - curr_hets - curr_homr;

		het_probs[mid] = 1.0;
		double sum = het_probs[mid];
		for (curr_hets = mid; curr_hets > 1; curr_hets -= 2){
			het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
			sum += het_probs[curr_hets - 2];

			// 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
			++curr_homr;
			++curr_homc;
		}

		curr_hets = mid;
		curr_homr = (rare_copies - mid) / 2;
		curr_homc = genotypes - curr_hets - curr_homr;
		for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2){
			het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /((curr_hets + 2.0) * (curr_hets + 1.0));
			sum += het_probs[curr_hets + 2];

			// add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
			--curr_homr;
			--curr_homc;
		}

		for (i = 0; i <= rare_copies; i++) het_probs[i] /= sum;

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
