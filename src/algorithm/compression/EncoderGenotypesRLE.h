#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>

#include "../../io/bcf/BCFEntry.h"
#include "../../core/base/EntryHotMeta.h"
#include "../../core/base/GTRecords.h"
#include "../../core/StreamContainer.h"

namespace Tachyon{
namespace Algorithm{

#define PACK_RLE_BIALLELIC(A, B, SHIFT, ADD) BCF::BCF_UNPACK_GENOTYPE(A) << (SHIFT + ADD) | BCF::BCF_UNPACK_GENOTYPE(B) << (ADD) | (A & ADD)
#define PACK_RLE_SIMPLE(A, B, SHIFT) ((A >> 1) << (SHIFT + 1)) | ((B >> 1) << 1) | (A & 1)

// Todo: These should all be moved to VCF or BCF entry class
struct EncoderGenotypesRLEHelper{
	typedef EncoderGenotypesRLEHelper self_type;

	EncoderGenotypesRLEHelper(void) :
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

	EncoderGenotypesRLEHelper(const U64 expectedSamples) :
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
	~EncoderGenotypesRLEHelper(){}

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
	// Modified to use Tomahawk data
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

class EncoderGenotypesRLE {
	typedef EncoderGenotypesRLE self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef BCF::BCFEntry bcf_type;
	typedef Core::EntryHotMetaBase meta_base_type;
	typedef EncoderGenotypesRLEHelper helper_type;
	typedef Core::StreamContainer container_type;

	typedef struct __RLEAssessHelper{
		explicit __RLEAssessHelper(void) : mixedPhasing(1), hasMissing(1), word_width(1), n_runs(0){}
		__RLEAssessHelper(const BYTE& word_width,
				          const U64& n_runs,
						  const bool& mixedPhasing,
						  const bool& hasMissing) :
			mixedPhasing(mixedPhasing),
			hasMissing(hasMissing),
			word_width(word_width),
			n_runs(n_runs)
		{}
		~__RLEAssessHelper(){}

		bool mixedPhasing;
		bool hasMissing;
		BYTE word_width;
		U64 n_runs;

	} rle_helper_type;

public:
	EncoderGenotypesRLE() :
		n_samples(0)
	{
	}

	EncoderGenotypesRLE(const U64 samples) :
		n_samples(samples),
		helper(samples)
	{
	}

	~EncoderGenotypesRLE(){}
	bool Encode(const bcf_type& line, meta_base_type& meta_base, container_type& runs, container_type& simple, U64& n_runs, const U32* const ppa);
	inline void setSamples(const U64 samples){ this->n_samples = samples; this->helper.setExpected(samples); }

private:
	const rle_helper_type assessRLEBiallelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessRLEnAllelic(const bcf_type& line, const U32* const ppa);
	template <class T> bool EncodeSimple(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeRLE(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa, const bool hasMissing = true, const bool hasMixedPhase = true);
	template <class T> bool EncodeRLESimple(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);

private:
	U64 n_samples;             // number of samples
	helper_type helper;        // support stucture
};

template <class T>
bool EncoderGenotypesRLE::EncodeSimple(const bcf_type& line, container_type& simple, U64& n_runs, const U32* const ppa){
	BYTE shift_size = 3;
	if(sizeof(T) == 2) shift_size = 7;
	if(sizeof(T) == 4) shift_size = 15;

	// Virtual byte offset into start of genotypes
	// in BCF entry
	U32 internal_pos = line.p_genotypes;

	// Pack genotypes as
	// allele A | alleleB | isPhased
	if(sizeof(T) <= 2){
		for(U32 i = 0; i < this->n_samples * 2; i += 2){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]+1]);

			//const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
			//const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
			const T packed = ((fmt_type_value2 >> 1) << (shift_size + 1)) |
							 ((fmt_type_value1 >> 1) << 1) |
							 (fmt_type_value2 & 1);
			simple += packed;
		}
	} else if(sizeof(T) == 4){
		for(U32 i = 0; i < this->n_samples * 4; i += 4){
			//const S16& fmt_type_value1 = *reinterpret_cast<const S16* const>(&line.data[internal_pos++]);
			//const S16& fmt_type_value2 = *reinterpret_cast<const S16* const>(&line.data[internal_pos++]);
			const S16& fmt_type_value1 = *reinterpret_cast<const S16* const>(&line.data[internal_pos+4*ppa[0]]);
			const S16& fmt_type_value2 = *reinterpret_cast<const S16* const>(&line.data[internal_pos+4*ppa[0]+2]);
			const T packed = ((fmt_type_value2 >> 1) << (shift_size + 1)) |
							 ((fmt_type_value1 >> 1) << 1) |
							 (fmt_type_value2 & 1);
			simple += packed;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal number of alleles..." << std::endl;
		exit(1);
	}
	n_runs = this->n_samples;
	simple.n_additions += n_runs;

	return(true);
}

template <class T>
bool EncoderGenotypesRLE::EncodeRLE(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa, const bool hasMissing, const bool hasMixedPhase){
	U32 internal_pos = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength = 0;
	T length = 1;
	T RLE = 0;
	// pack
	// UNPACK(value) << shift +1 | UNPACK(value2) << 1 | phase
	// shift = 2 if hasMissing == TRUE
	// shift = 1 if HasMissing == FALSE
	// add = 1 if hasMixedPhase == TRUE
	const BYTE shift = hasMissing ? 2 : 1;
	const BYTE add = hasMixedPhase ? 1 : 0;

	// Run limits
	const T limit = pow(2, 8*sizeof(T) - (2*(1+hasMissing)+hasMixedPhase)) - 1;

	// Genotype maps
	// Map to 0,1,4,5
	const BYTE*    map = Constants::ALLELE_REDUCED_MAP;
	if(shift == 2) map = Constants::ALLELE_SELF_MAP;

	//std::cerr << ppa[0] << std::endl;
	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]+1]);
	T packed = PACK_RLE_BIALLELIC(fmt_type_value2, fmt_type_value1, shift, add);
	//std::cout << (U32)packed << '\t';

	// MSB contains phasing information
	this->helper.phased = (fmt_type_value2 & 1);

	U32 j = 1;
	BYTE last_phase = (fmt_type_value2 & 1);
	for(U32 i = 2; i < this->n_samples * 2; i += 2, ++j){
		//std::cerr << ppa[j] << std::endl;
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[j]]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[j]+1]);
		const T packed_internal = PACK_RLE_BIALLELIC(fmt_type_value2, fmt_type_value1, shift, add);
		last_phase = (fmt_type_value2 & 1);

		//std::cout << (U32)packed_internal << '\t';

		if(packed != packed_internal || length == limit){
			// Prepare RLE
			RLE = length;
			RLE <<= (2*shift + add);
			if(add) RLE |= (fmt_type_value2 & 1);
			assert((RLE >> (2*shift + add)) == length);

			// Set meta phased flag bit
			if((fmt_type_value2 & 1) != 1) this->helper.phased = false;

			// Push RLE to buffer
			runs += RLE;

			// Update genotype and allele counts
			this->helper[map[packed >> add]] += length;
			this->helper.countsAlleles[packed >> (shift + add)] += length;
			this->helper.countsAlleles[(packed >> add) & ((1 << shift) - 1)]  += length;

			// Reset and update
			sumLength += length;
			length = 0;
			packed = packed_internal;
			++n_runs;
		}
		++length;
	}
	// Last entry
	// Prepare RLE
	RLE = length;
	RLE <<= 2 * shift + add;
	if(add) RLE |= (last_phase & 1);
	assert((RLE >> (2*shift + add)) == length);

	// Set meta phased flag bit
	if((packed & 1) != 1) this->helper.phased = false;

	// Push RLE to buffer
	runs += RLE;
	++n_runs;

	// Update genotype and allele counts
	this->helper[map[packed >> add]] += length;
	this->helper.countsAlleles[packed >> (shift + add)] += length;
	this->helper.countsAlleles[(packed >> add) & ((1 << shift) - 1)]  += length;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	runs.n_additions += n_runs;

	// Calculate basic stats
	this->helper.calculateMGF();
	this->helper.calculateHardyWeinberg();
	return(true);
}

template <class T>
bool EncoderGenotypesRLE::EncodeRLESimple(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa){
	U32 internal_pos = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength = 0;
	T length = 1;
	T RLE = 0;

	const BYTE shift = ceil(log2(line.body->n_allele + 1));

	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]+1]);

	//const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
	//const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
	T packed = PACK_RLE_SIMPLE(fmt_type_value2, fmt_type_value1, shift);

	// MSB contains phasing information
	this->helper.phased = (fmt_type_value2 & 1);

	// Run limits
	const T limit = pow(2, 8*sizeof(T) - (2*shift+1)) - 1;

	for(U32 i = 2; i < this->n_samples * 2; i += 2){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos+2*ppa[0]+1]);

		//const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
		//const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos++]);
		const T packed_internal = PACK_RLE_SIMPLE(fmt_type_value2, fmt_type_value1, shift);

		if(packed != packed_internal || length == limit){
			// Prepare RLE
			RLE = length;
			RLE <<= 2*shift + 1;
			RLE |= packed;

			//if(sizeof(T) == 1) std::cerr << line.body->POS+1 << ": " << std::bitset<sizeof(T)*8>(RLE) << '\t' << std::bitset<sizeof(T)*8>(packed) << std::endl;

			// Set meta phased flag bit
			if((fmt_type_value2 & 1) != 1) this->helper.phased = false;
			assert((RLE >> (2*shift + 1)) == length);

			// Push RLE to buffer
			runs += RLE;

			// Reset and update
			sumLength += length;
			length = 0;
			packed = packed_internal;
			++n_runs;
		}
		++length;
	}
	// Last entry
	// Prepare RLE
	RLE = length;
	RLE <<= 2*shift + 1;
	RLE |= packed;
	assert((RLE >> (2*shift + 1)) == length);
	//if(sizeof(T) == 1){ std::cerr << line.body->POS+1 << ": " << std::bitset<sizeof(T)*8>(RLE) << std::endl; exit(1);}

	// Set meta phased flag bit
	if((packed & 1) != 1) this->helper.phased = false;

	// Push RLE to buffer
	runs += RLE;
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);

	runs.n_additions += n_runs;

	return(true);
}

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
