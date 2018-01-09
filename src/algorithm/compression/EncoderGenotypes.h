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

#define PACK_RLE_BIALLELIC(A, B, SHIFT, ADD) BCF::BCF_UNPACK_GENOTYPE(A) << (SHIFT + ADD) | BCF::BCF_UNPACK_GENOTYPE(B) << (ADD) | (A & ADD)
#define PACK_RLE_SIMPLE(A, B, SHIFT) ((A >> 1) << (SHIFT + 1)) | ((B >> 1) << 1) | (A & 1)

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
				mixedPhasing(1),
				hasMissing(1),
				word_width(1),
				n_runs(0)
		{}
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
	EncoderGenotypes();
	EncoderGenotypes(const U64 samples);
	~EncoderGenotypes();
	bool Encode(const bcf_type& line, meta_type& meta_base, container_type& runs, container_type& simple, container_type& support, U64& n_runs, const U32* const ppa);
	inline void setSamples(const U64 samples){ this->n_samples = samples; }

private:
	const rle_helper_type assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessMploidRLEBiallelic(const bcf_type& line, const U32* const ppa);
	const rle_helper_type assessMploidRLEnAllelic(const bcf_type& line, const U32* const ppa);

	template <class T> bool EncodeDiploidSimple(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeDiploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa, const bool hasMissing = true, const bool hasMixedPhase = true);
	template <class T> bool EncodeDiploidRLEnAllelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeMploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);
	template <class T> bool EncodeMploidRLENallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa);

private:
	U64 n_samples; // number of samples
};

template <class T>
bool EncoderGenotypes::EncodeDiploidSimple(const bcf_type& line, container_type& simple, U64& n_runs, const U32* const ppa){
	BYTE shift_size = 3;
	if(sizeof(T) == 2) shift_size = 7;
	if(sizeof(T) == 4) shift_size = 15;

	// Virtual byte offset into start of genotypes
	// in BCF entry
	U32 internal_pos = line.p_genotypes;

	// Pack genotypes as
	// allele A | alleleB | isPhased
	if(sizeof(T) <= 2){
		U32 j = 0;
		for(U32 i = 0; i < this->n_samples * 2; i += 2, ++j){
			const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]]);
			const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]+sizeof(SBYTE)]);

			const T packed = ((allele2 >> 1) << (shift_size + 1)) |
							 ((allele1 >> 1) << 1) |
							 (allele2 & 1);
			simple += packed;
		}
	} else if(sizeof(T) == 4){
		U32 j = 0;
		for(U32 i = 0; i < this->n_samples * 4; i += 4, ++j){
			const S16& allele1 = *reinterpret_cast<const S16* const>(&line.data[internal_pos + 2*sizeof(S16)*ppa[j]]);
			const S16& allele2 = *reinterpret_cast<const S16* const>(&line.data[internal_pos + 2*sizeof(S16)*ppa[j]+sizeof(S16)]);
			const T packed = ((allele2 >> 1) << (shift_size + 1)) |
							 ((allele1 >> 1) << 1) |
							 (allele2 & 1);
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
bool EncoderGenotypes::EncodeDiploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa, const bool hasMissing, const bool hasMixedPhase){
	U32 internal_pos = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength = 0;
	T length = 1;
	T RLE = 0;
	// pack
	// UNPACK(value) << shift +1 | UNPACK(value2) << 1 | phase
	// shift = 2 if hasMissing == TRUE
	// shift = 1 if HasMissing == FALSE
	// add = 1 if hasMixedPhase == TRUE
	const BYTE shift = hasMissing    ? 2 : 1;
	const BYTE add   = hasMixedPhase ? 1 : 0;

	// Run limits
	const T limit = pow(2, 8*sizeof(T) - (2*(1 + hasMissing) + hasMixedPhase)) - 1;

	// Genotype maps
	// Map to 0,1,4,5
	const BYTE*    map = Constants::ALLELE_REDUCED_MAP;
	if(shift == 2) map = Constants::ALLELE_SELF_MAP;

	//std::cerr << ppa[0] << std::endl;
	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[0]+sizeof(SBYTE)]);
	T packed = PACK_RLE_BIALLELIC(allele2, allele1, shift, add);
	//std::cout << (U32)packed << '\t';

	// MSB contains phasing information
	//this->helper.phased = (fmt_type_value2 & 1);

	U32 j = 1;
	BYTE last_phase = (allele2 & 1);
	for(U32 i = 2; i < this->n_samples * 2; i += 2, ++j){
		//std::cerr << ppa[j] << std::endl;
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]+sizeof(SBYTE)]);
		const T packed_internal = PACK_RLE_BIALLELIC(allele2, allele1, shift, add);
		last_phase = (allele2 & 1);

		//std::cout << (U32)packed_internal << '\t';

		if(packed != packed_internal || length == limit){
			// Prepare RLE
			RLE = length;
			RLE <<= (2*shift + add);
			if(add) RLE |= (allele2 & 1);
			assert((RLE >> (2*shift + add)) == length);

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
	RLE <<= 2 * shift + add;
	if(add) RLE |= (last_phase & 1);
	assert((RLE >> (2*shift + add)) == length);

	// Push RLE to buffer
	runs += RLE;
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	runs.n_additions += n_runs;

	return(true);
}

template <class T>
bool EncoderGenotypes::EncodeDiploidRLEnAllelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa){
	U32 internal_pos = line.p_genotypes; // virtual byte offset of genotype start
	U32 sumLength = 0;
	T length = 1;
	T RLE = 0;

	const BYTE shift = ceil(log2(line.body->n_allele + 1));

	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[0]+sizeof(SBYTE)]);
	T packed = PACK_RLE_SIMPLE(allele2, allele1, shift);

	// Run limits
	const T limit = pow(2, 8*sizeof(T) - (2*shift+1)) - 1;

	U32 j = 1;
	for(U32 i = 2; i < this->n_samples * 2; i += 2, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos + 2*sizeof(SBYTE)*ppa[j]+sizeof(SBYTE)]);
		const T packed_internal = PACK_RLE_SIMPLE(allele2, allele1, shift);

		if(packed != packed_internal || length == limit){
			// Prepare RLE
			RLE = length;
			RLE <<= 2*shift + 1;
			RLE |= packed;

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
