#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>

#include "../../containers/variantblock.h"
#include "../../core/meta_hot.h"
#include "../../io/bcf/BCFEntry.h"
#include "../../core/genotype_summary.h"

namespace tachyon{
namespace algorithm{

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

#define YON_PACK_GT_DIPLOID(A, B, SHIFT, ADD) (bcf::BCF_UNPACK_GENOTYPE(A) << ((SHIFT) + (ADD))) | (bcf::BCF_UNPACK_GENOTYPE(B) << (ADD)) | ((A) & (ADD))
#define YON_PACK_GT_DIPLOID_NALLELIC(A, B, SHIFT, ADD) (((A) >> 1) << ((SHIFT) + (ADD))) | (((B) >> 1) << (ADD)) | ((A) & (ADD))

class GenotypeEncoder {
private:
	typedef GenotypeEncoder              self_type;
	typedef io::BasicBuffer              buffer_type;
	typedef bcf::BCFEntry                bcf_type;
	typedef core::MetaEntry              meta_type;
	typedef containers::GenotypesSummary helper_type;
	typedef containers::DataContainer    container_type;
	typedef containers::VariantBlock     block_type;

	typedef struct __RLEAssessHelper{
		explicit __RLEAssessHelper(void) :
				word_width(1),
				n_runs(0)
		{}
		__RLEAssessHelper(const BYTE& word_width,
				          const U64& n_runs) :
			word_width(word_width),
			n_runs(n_runs)
		{}
		~__RLEAssessHelper(){}

		BYTE word_width;
		U64 n_runs;

	} rle_helper_type;

public:
	GenotypeEncoder();
	GenotypeEncoder(const U64 samples);
	~GenotypeEncoder();
	bool Encode(const bcf_type& line, meta_type& meta_base, block_type& block, const U32* const ppa);
	inline void setSamples(const U64 samples){ this->n_samples = samples; }

private:
	const rle_helper_type assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa) const;
	const rle_helper_type assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa) const;
	const rle_helper_type assessMploidRLEBiallelic(const bcf_type& line, const U32* const ppa) const;
	const rle_helper_type assessMploidRLEnAllelic(const bcf_type& line, const U32* const ppa) const;

	template <class YON_STORE_TYPE, class BCF_GT_TYPE = SBYTE> bool EncodeBCFStyle(const bcf_type& line, container_type& container, U64& n_runs) const;
	template <class YON_RLE_TYPE, class BCF_GT_TYPE = SBYTE> bool EncodeDiploidBCF(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa) const;
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEBiallelic(const bcf_type& line, container_type& runs, const U32* const ppa, const rle_helper_type& helper) const;
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEnAllelic(const bcf_type& line, container_type& runs, const U32* const ppa, const rle_helper_type& helper) const;
	template <class T> bool EncodeMploidRLEBiallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa) const;
	template <class T> bool EncodeMploidRLENallelic(const bcf_type& line, container_type& runs, U64& n_runs, const U32* const ppa) const;

private:
	U64 n_samples; // number of samples
};

template <class YON_STORE_TYPE, class BCF_GT_TYPE>
bool GenotypeEncoder::EncodeBCFStyle(const bcf_type& line,
                                     container_type& simple,
                                                U64& n_runs) const
{
	const BYTE ploidy = line.gt_support.ploidy;
	U32 bcf_gt_pos = line.formatID[0].l_offset;
	const BCF_GT_TYPE missing_value = (BCF_GT_TYPE)1 << (sizeof(BCF_GT_TYPE)*8 - 1);
	const BCF_GT_TYPE EOV_value     = missing_value + 1;


	// Pack genotypes as
	// allele | phasing
	U32 j = 0;
	for(U32 i = 0; i < this->n_samples * ploidy; i += ploidy, ++j){
		for(U32 p = 0; p < ploidy; ++p){
			const BCF_GT_TYPE& allele = *reinterpret_cast<const BCF_GT_TYPE* const>(&line.data[bcf_gt_pos]);
			if((allele >> 1) == 0) simple.AddLiteral((YON_STORE_TYPE)0); // missing
			else if(allele == EOV_value ) simple.AddLiteral((YON_STORE_TYPE)1); // eov
			else { // otherwise
				// Add 1 because 1 is reserved for EOV
				const YON_STORE_TYPE val = ((allele >> 1) + 1) << 1 | (allele & 1);
				simple.AddLiteral(val);
			}
			bcf_gt_pos += sizeof(BCF_GT_TYPE);
		}
	}

	n_runs = this->n_samples*ploidy;
	simple.n_additions += n_runs;

	return(true);
}

template <class YON_RLE_TYPE, class BCF_GT_TYPE>
bool GenotypeEncoder::EncodeDiploidBCF(const bcf_type& line,
		                               container_type& simple,
										          U64& n_runs,
                                     const U32* const  ppa) const
{
	const BYTE ploidy = 2;
	BYTE shift_size = 3;
	if(sizeof(YON_RLE_TYPE) == 2) shift_size = 7;
	if(sizeof(YON_RLE_TYPE) == 4) shift_size = 15;

	// Virtual byte offset into start of genotypes
	// in BCF entry
	U32 bcf_gt_pos = line.formatID[0].l_offset;

	// Pack genotypes as
	// allele A | alleleB | isPhased
	U32 j = 0;
	for(U32 i = 0; i < this->n_samples * ploidy; i += ploidy, ++j){
		const BCF_GT_TYPE& allele1 = *reinterpret_cast<const BCF_GT_TYPE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BCF_GT_TYPE)*ppa[j]]);
		const BCF_GT_TYPE& allele2 = *reinterpret_cast<const BCF_GT_TYPE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BCF_GT_TYPE)*ppa[j]+sizeof(BCF_GT_TYPE)]);

		const YON_RLE_TYPE packed = ((allele2 >> 1) << (shift_size + 1)) |
				                    ((allele1 >> 1) << 1) |
									 (allele2 &  1);
		simple.AddLiteral(packed);
	}

	n_runs = this->n_samples;
	simple.n_additions += n_runs;

	return(true);
}

template <class YON_RLE_TYPE>
bool GenotypeEncoder::EncodeDiploidRLEBiallelic(const bcf_type& line,
		                                         container_type& runs,
												const U32* const ppa,
								          const rle_helper_type& helper) const
{
	const BYTE ploidy   = 2;
	U32 bcf_gt_pos      = line.formatID[0].l_offset; // virtual byte offset of genotype start
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	const BYTE shift    = line.gt_support.hasMissing    ? 2 : 1;
	const BYTE add      = line.gt_support.mixedPhasing  ? 1 : 0;

	// temp
	//const U64 runs_start_pos = runs.buffer_data_uncompressed.pointer;

	// Run limits
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (ploidy*(1 + line.gt_support.hasMissing) + line.gt_support.mixedPhasing)) - 1;

	// Genotype maps
	// Map to 0,1,4,5
	//const BYTE*    map = Constants::ALLELE_REDUCED_MAP;
	//if(shift == 2) map = Constants::ALLELE_SELF_MAP;

	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[0]+sizeof(SBYTE)]);
	YON_RLE_TYPE packed = YON_PACK_GT_DIPLOID(allele2, allele1, shift, add);

	U32 ppa_pos = 1;
	BYTE last_phase = (allele2 & 1);
	U64 n_runs = 0;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++ppa_pos){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[ppa_pos]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(SBYTE)*ppa[ppa_pos]+sizeof(SBYTE)]);
		const YON_RLE_TYPE packed_internal = YON_PACK_GT_DIPLOID(allele2, allele1, shift, add);
		last_phase = (allele2 & 1);

		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = length;
			RLE <<= (ploidy*shift + add);
			RLE |= packed;
			assert((RLE >> (ploidy*shift + add)) == length);

			// Push RLE to buffer
			runs.AddLiteral(RLE);

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
	runs.AddLiteral(RLE);
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	assert(helper.n_runs == n_runs);
	runs.n_additions += n_runs;

#if ENCODER_GT_DEBUG == 1
	std::cout << 0 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif
	return(true);
}

template <class YON_RLE_TYPE>
bool GenotypeEncoder::EncodeDiploidRLEnAllelic(const bcf_type& line,
		                                       container_type& runs,
											 const U32* const  ppa,
										const rle_helper_type& helper) const
{
	const BYTE ploidy   = 2;
	U32 bcf_gt_pos      = line.formatID[0].l_offset; // virtual byte offset of genotype start
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	const BYTE shift    = ceil(log2(line.body->n_allele + line.gt_support.hasMissing + line.gt_support.hasEOV + 1)); // Bits occupied per allele
	const BYTE add      = line.gt_support.mixedPhasing  ? 1 : 0;
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE)  - (ploidy*shift + line.gt_support.mixedPhasing)) - 1;

	// Setup first run
	BYTE allele1 = *reinterpret_cast<const BYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BYTE)*ppa[0]]);
	BYTE allele2 = *reinterpret_cast<const BYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BYTE)*ppa[0]+sizeof(SBYTE)]);
	if(allele1 == 0x81){
		allele1 = 1;
		//std::cerr << "eov" << std::endl;
	}
	if(allele2 == 0x81){
		allele2 = 1;
		//std::cerr << "eov" << std::endl;
	}

	YON_RLE_TYPE packed  = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add);

	U32 j = 1;
	U64 n_runs = 0;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		BYTE allele1 = *reinterpret_cast<const BYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BYTE)*ppa[j]]);
		BYTE allele2 = *reinterpret_cast<const BYTE* const>(&line.data[bcf_gt_pos + ploidy*sizeof(BYTE)*ppa[j]+sizeof(SBYTE)]);
		if(allele1 == 0x81){
			allele1 = 1;
			//std::cerr << "eov" << std::endl;
		}
		if(allele2 == 0x81){
			allele2 = 1;
			//std::cerr << "eov" << std::endl;
		}

		const YON_RLE_TYPE packed_internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add);


		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = (length << (ploidy*shift + add)) | packed;

			//std::cerr << (int)line.body->n_allele << '\t' << std::bitset<32>(RLE) << '\t' << std::bitset<32>(RLE >> (ploidy*shift + add)) << '\t' << std::bitset<32>(length) << std::endl;
			assert((RLE >> (ploidy*shift + add)) == length);
			assert(length != 0);

			// Push RLE to buffer
			runs.AddLiteral(RLE);

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
	RLE = (length << (ploidy*shift + add)) | packed;
	//std::cerr << "last:" << (int)line.body->n_allele << '\t' << std::bitset<32>(RLE) << '\t' << std::bitset<32>(RLE >> (ploidy*shift + add)) << '\t' << std::bitset<32>(length) << std::endl;
	assert((RLE >> (ploidy*shift + add)) == length);
	assert(length != 0);

	// Push RLE to buffer
	runs.AddLiteral(RLE);
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	assert(helper.n_runs == n_runs);

	runs.n_additions += n_runs;
#if ENCODER_GT_DEBUG == 1
	std::cout << 1 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif

	return(true);
}

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
