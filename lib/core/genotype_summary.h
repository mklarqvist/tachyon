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

struct GenotypeSummaryObject{
	GenotypeSummaryObject() : counts_(0), countsA_(0), countsB_(0){}
	~GenotypeSummaryObject(){}

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

struct GenotypeSummary{
	typedef U64             value_type;
	typedef core::MetaEntry meta_type;

public:
	GenotypeSummary() :
		n_alleles_(10),
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

	GenotypeSummary(const BYTE n_alleles) :
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

	~GenotypeSummary(){
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

	U64 alleleCountA(void) const{
		U64 total = 0;
		for(U32 i = 2; i < this->n_alleles_; ++i){
			total += this->vectorA_[i];
		}
		return(total);
	}

	U64 alleleCountB(void) const{
		U64 total = 0;
		for(U32 i = 2; i < this->n_alleles_; ++i){
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

	/**<
	// This code implements an exact SNP test of Hardy-Weinberg Equilibrium as described in
	// Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
	// Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
	//
	// Written by Jan Wigginton
	// Modified to use Tachyon data by Marcus D. R. Klarqvist (https;//github.com/mklarqvist/tachyon)
	*/
	std::vector<double> calculateHardyWeinberg(const meta_type& meta) const{
		if(meta.n_alleles == 1) return std::vector<double>();
		const BYTE n_limit = meta.n_alleles + 2 > this->n_alleles_ ? this->n_alleles_ - 1 : meta.n_alleles - 1;
		std::vector<double> results(n_limit, 1);

		for(U32 i = 0; i < n_limit; ++i)
			results[i] = this->__calculateHardyWeinberg(2, 3+i);

		return(results);
	}

	std::vector<double> calculateAlleleFrequency(const meta_type& meta) const{
		const BYTE n_limit = meta.n_alleles + 2 > this->n_alleles_
		                   ? this->n_alleles_ - 2
		                   : meta.n_alleles;
		std::vector<double> results(n_limit, 0);

		U64 n_total = 0;
		for(U32 i = 0; i < n_limit; ++i){
			results[i] = this->vectorA_[2+i] + this->vectorB_[2+i];
			n_total   += results[i];
		}
		for(U32 i = 0; i < n_limit; ++i) results[i] /= n_total;

		return(results);
	}

	std::vector<double> calculateInbreedingCoefficient(const meta_type& meta) const{
		// Allele frequency of A
		const double p = ((double)2*this->matrix_[2][2] + this->matrix_[2][3] + this->matrix_[3][2]) / (2*this->genotypeCount());
		// Genotype frequency of heterozygotes
		const double pg = ((double)this->matrix_[2][3] + this->matrix_[3][2]) / this->genotypeCount();
		// Expected heterozygosity
		const double exp = 2*p*(1-p);
		// Population inbreeding coefficient: F
		const double f_pic = exp > 0 ? (exp-pg)/exp : 0;

		return(std::vector<double>());
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
		const BYTE shift    = ceil(log2(gt_simple_container.getMeta().getNumberAlleles() + 2 + 1)); // Bits occupied per allele, 1 value for missing
		const BYTE add      = gt_simple_container.getMeta().isGTMixedPhasing() ? 1 : 0;
		//const BYTE matrix_add = !gt_simple_container.getMeta().isMixedPloidy();

		if(gt_simple_container.getMeta().n_alleles + 2 > this->n_alleles_){
			std::cerr << "too many alleles: " << gt_simple_container.getMeta().n_alleles + 2 << "/" << (int)this->n_alleles_ << std::endl;
			return;
		}

		for(U32 i = 0; i < gt_simple_container.size(); ++i){
			const U64 length   = YON_GT_RLE_LENGTH(gt_simple_container.at(i), shift, add);
			const BYTE alleleA = YON_GT_RLE_ALLELE_A(gt_simple_container.at(i), shift, add);
			const BYTE alleleB = YON_GT_RLE_ALLELE_B(gt_simple_container.at(i), shift, add);
			this->matrix_[alleleA][alleleB] += length;
			this->vectorA_[alleleA] += length;
			this->vectorB_[alleleB] += length;
		}
	}

	template <class T>
	inline void operator+=(const GenotypeContainerDiploidBCF<T>& gt_diploid_bcf_container){
		const BYTE shift = (sizeof(T)*8 - 1) / 2;

		for(U32 i = 0; i < gt_diploid_bcf_container.size(); ++i){
			const U16 alleleA = YON_GT_DIPLOID_BCF_A(gt_diploid_bcf_container.at(i), shift);
			const U16 alleleB = YON_GT_DIPLOID_BCF_B(gt_diploid_bcf_container.at(i), shift);

			this->matrix_[alleleA][alleleB] += 1;
			this->vectorA_[alleleA] += 1;
			this->vectorB_[alleleB] += 1;
		}
	}

	inline void operator+=(const bcf::BCFEntry& entry){
		if(entry.hasGenotypes == false) return;

		if(entry.body->n_allele + 2 > this->n_alleles_){
			std::cerr << "too many alleles: " << entry.body->n_allele + 2 << "/" << (int)this->n_alleles_ << std::endl;
			return;
		}

		U32 internal_pos = entry.formatID[0].l_offset;
		U32 k = 0;
		for(U32 i = 0; i < 2*entry.body->n_sample; i += 2, ++k){
			const BYTE& fmt_type_value1 = *reinterpret_cast<const BYTE* const>(&entry.data[internal_pos++]);
			const BYTE& fmt_type_value2 = *reinterpret_cast<const BYTE* const>(&entry.data[internal_pos++]);
			BYTE alleleA = fmt_type_value2 >> 1;
			BYTE alleleB = fmt_type_value1 >> 1;
			alleleA += (alleleA != 0 ? 1 : 0);
			alleleB += (alleleB != 0 ? 1 : 0);
			//std::cerr << (int)alleleA << "," << (int)alleleB << std::endl;
			if(!(alleleA < this->n_alleles_ && alleleB < this->n_alleles_)){
				std::cerr << entry.body->n_allele << std::endl;
				std::cerr << internal_pos << "/" << entry.l_data << std::endl;
				std::cerr << "pos: " << i << "/" << 2*entry.body->n_sample << "@" << k << std::endl;
				std::cerr << (int)alleleA << "," << (int)alleleB << std::endl;
				std::cerr << std::bitset<8>(fmt_type_value2) << "," << std::bitset<8>(fmt_type_value1) << std::endl;
				exit(1);
			}

			++this->matrix_[alleleA][alleleB];
			//++this->vectorA_[alleleA];
			//++this->vectorB_[alleleB];
		}
	}

private:
	double __calculateHardyWeinberg(const U32 ref_target, const U32 alt_target) const{
		U64 obs_hets = this->matrix_[ref_target][alt_target] + this->matrix_[alt_target][ref_target];
		U64 obs_hom1 = this->matrix_[ref_target][ref_target];
		U64 obs_hom2 = this->matrix_[alt_target][alt_target];
		if(obs_hets + obs_hom1 + obs_hom2 == 0) return 1;

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

		return(p_hwe);
	}

public:
	BYTE  n_alleles_; // number of alleles (including EOV and missing)
	value_type** matrix_;
	value_type*  vectorA_;
	value_type*  vectorB_;
};

}
}



#endif /* ENCODERGENOTYPESRLEHELPER_H_ */
