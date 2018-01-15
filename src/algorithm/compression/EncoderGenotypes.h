#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>

#include "../../io/bcf/BCFEntry.h"
#include "../../core/base/GTRecords.h"
#include "../../core/base/MetaHot.h"
#include "../../core/StreamContainer.h"
#include "EncoderGenotypesHelper.h"

namespace Tachyon{
namespace Encoding{

#define ENCODER_GT_DEBUG 0

/*
 RLE biallelic diploid
 +---+---+---+---+
 |  RLE  | A/B|P |
 +---+---+---+---+
 ^       ^
 |       |
 LENGTH  BIALLELIC GT

*/

/*
 Simple biallelic
 +---+---+---+---+---+---+
 |  A    |   B   |   P   |
 +---+---+---+---+---+---+
 ^       ^       ^
 |       |       |
 A       B       P

*/

/*

 sizeof(T1) + sizeof(T2)

 RLE simple biallelic
 +---+---+---+---+---+---+---+---+
 |  RLE  |   A   |   B   |   P   |
 +---+---+---+---+---+---+---+---+
 ^       ^       ^       ^
 |       |       |       |
 LENGTH  A       B       P

*/

#define PACK_DIPLOID_BIALLELIC(A, B, SHIFT, ADD) BCF::BCF_UNPACK_GENOTYPE(A) << (SHIFT + ADD) | BCF::BCF_UNPACK_GENOTYPE(B) << (ADD) | (A & ADD)
#define PACK_DIPLOID_NALLELIC(A, B, SHIFT) ((A >> 1) << (SHIFT + 1)) | ((B >> 1) << 1) | (A & 1)

class EncoderGenotypes {
private:
	typedef EncoderGenotypes self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef BCF::BCFEntry bcf_type;
	typedef Core::MetaHot meta_type;
	typedef EncoderGenotypesHelper helper_type;
	typedef Core::StreamContainer container_type;

	typedef struct __RLEAssessHelper{
		explicit __RLEAssessHelper(void) :
				phased(0),
				mixedPhasing(1),
				hasMissing(1),
				anyNA(0),
				word_width(1),
				n_runs(0)
		{}
		__RLEAssessHelper(const BYTE& word_width,
				          const U64& n_runs,
						  const bool& phase,
						  const bool& mixedPhasing,
						  const bool& hasMissing,
						  const bool& anyNA) :
			phased(phase),
			mixedPhasing(mixedPhasing),
			hasMissing(hasMissing),
			anyNA(anyNA),
			word_width(word_width),
			n_runs(n_runs)
		{}
		~__RLEAssessHelper(){}

		bool phased;
		bool mixedPhasing;
		bool hasMissing;
		bool anyNA;
		BYTE word_width;
		U64 n_runs;

	} rle_helper_type;

public:
	EncoderGenotypes();
	EncoderGenotypes(const U64 samples);
	~EncoderGenotypes();
	bool Encode(const bcf_type& line, meta_type& meta_base, container_type& runs, container_type& simple, container_type& support, const U32* const ppa);
	inline void setSamples(const U64 samples){ this->n_samples = samples; }

private:
	const rle_helper_type assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessMploidRLEBiallelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessMploidRLEnAllelic(const bcf_type& line, const U32* const ppa);

	template <class YON_RLE_TYPE, class BCF_GT_TYPE = SBYTE> bool EncodeDiploidBCF(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa, const bool hasMissing = true, const bool hasMixedPhase = true);
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEnAllelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeMploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeMploidRLENallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);

private:
	U64 n_samples; // number of samples
};

template <class YON_RLE_TYPE, class BCF_GT_TYPE>
bool EncoderGenotypes::EncodeDiploidBCF(const bcf_type& line,
		                                   container_type& simple,
										              U64& n_runs,
										  const U32* const ppa){
	const BYTE ploidy = 2;
	BYTE shift_size = 3;
	if(sizeof(YON_RLE_TYPE) == 2) shift_size = 7;
	if(sizeof(YON_RLE_TYPE) == 4) shift_size = 15;

	// Virtual byte offset into start of genotypes
	// in BCF entry
	U32 bcf_gt_pos = line.p_genotypes;

	// Pack genotypes as
	// allele A | alleleB | isPhased
	U32 j = 0;
	for(U32 i = 0; i < this->n_samples * ploidy; i += ploidy, ++j){
		const BCF_GT_TYPE& allele1 = *reinterpret_cast<const BCF_GT_TYPE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BCF_GT_TYPE)*ppa[j]]);
		const BCF_GT_TYPE& allele2 = *reinterpret_cast<const BCF_GT_TYPE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BCF_GT_TYPE)*ppa[j]+sizeof(BCF_GT_TYPE)]);

		const YON_RLE_TYPE packed = ((allele2 >> 1) << (shift_size + 1)) |
				                    ((allele1 >> 1) << 1) |
									(allele2 & 1);
		simple += packed;
	}

	n_runs = this->n_samples;
	simple.n_additions += n_runs;

	return(true);
}

template <class YON_RLE_TYPE>
bool EncoderGenotypes::EncodeDiploidRLEBiallelic(const bcf_type& line,
		                                         container_type& runs,
												            U64& n_runs,
												const U32* const ppa,
												      const bool hasMissing,
													  const bool hasMixedPhase){
	const BYTE ploidy   = 2;
	U32 bcf_gt_pos      = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	// pack
	// UNPACK(value) << shift +1 | UNPACK(value2) << 1 | phase
	// shift = 2 if hasMissing == TRUE
	// shift = 1 if HasMissing == FALSE
	// add = 1 if hasMixedPhase == TRUE
	const BYTE shift    = hasMissing    ? 2 : 1;
	const BYTE add      = hasMixedPhase ? 1 : 0;

	// temp
	//const U64 runs_start_pos = runs.buffer_data_uncompressed.pointer;

	// Run limits
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (ploidy*(1 + hasMissing) + hasMixedPhase)) - 1;

	// Genotype maps
	// Map to 0,1,4,5
	//const BYTE*    map = Constants::ALLELE_REDUCED_MAP;
	//if(shift == 2) map = Constants::ALLELE_SELF_MAP;

	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]+sizeof(SBYTE)]);
	YON_RLE_TYPE packed = PACK_DIPLOID_BIALLELIC(allele2, allele1, shift, add);

	U32 ppa_pos = 1;
	BYTE last_phase = (allele2 & 1);
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++ppa_pos){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[ppa_pos]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[ppa_pos]+sizeof(SBYTE)]);
		const YON_RLE_TYPE packed_internal = PACK_DIPLOID_BIALLELIC(allele2, allele1, shift, add);
		last_phase = (allele2 & 1);

		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = length;
			RLE <<= (ploidy*shift + add);
			RLE |= packed;
			assert((RLE >> (ploidy*shift + add)) == length);

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
	RLE <<= ploidy * shift + add;
	if(add) RLE |= (last_phase & 1);
	assert((RLE >> (ploidy*shift + add)) == length);

	// Push RLE to buffer
	runs += RLE;
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	runs.n_additions += n_runs;

	/*
	if(add == 0){
		Core::TachyonRunNoPhase<YON_RLE_TYPE, false>* c = reinterpret_cast< Core::TachyonRunNoPhase<YON_RLE_TYPE, false>* >(&runs.buffer_data_uncompressed[runs_start_pos]);
		for(U32 i = 0; i < n_runs; ++i){
			std::cout << (int)c->alleleA << '|' << (int)c->alleleB << ':' << (int)c->runs << '\t';
			++c;
		}
		std::cout << std::endl;
	}
	*/


#if ENCODER_GT_DEBUG == 1
	std::cout << 0 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif
	return(true);
}

template <class YON_RLE_TYPE>
bool EncoderGenotypes::EncodeDiploidRLEnAllelic(const bcf_type& line,
		                                        container_type& runs,
												           U64& n_runs,
													 const U32* const ppa){
	const BYTE ploidy   = 2;
	U32 bcf_gt_pos      = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	const BYTE shift    = ceil(log2(line.body->n_allele + 1)); // Bits occupied per allele
	// Run limits
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (ploidy*shift + 1)) - 1;

	// First
	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]+sizeof(SBYTE)]);
	YON_RLE_TYPE packed = PACK_DIPLOID_NALLELIC(allele2, allele1, shift);

	U32 j = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[j]+sizeof(SBYTE)]);
		const YON_RLE_TYPE packed_internal = PACK_DIPLOID_NALLELIC(allele2, allele1, shift);

		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = length;
			RLE <<= ploidy*shift + 1;
			RLE |= packed;

			assert((RLE >> (ploidy*shift + 1)) == length);

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
	RLE <<= ploidy*shift + 1;
	RLE |= packed;
	assert((RLE >> (ploidy*shift + 1)) == length);

	// Push RLE to buffer
	runs += RLE;
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);

	runs.n_additions += n_runs;
#if ENCODER_GT_DEBUG == 1
	std::cout << 1 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif

	return(true);
}

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
