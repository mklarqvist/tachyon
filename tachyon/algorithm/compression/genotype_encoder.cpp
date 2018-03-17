#include "genotype_encoder.h"

namespace tachyon{
namespace algorithm{

GenotypeEncoder::GenotypeEncoder() :
	n_samples(0)
{
}

GenotypeEncoder::GenotypeEncoder(const U64 samples) :
	n_samples(samples)
{
}

GenotypeEncoder::~GenotypeEncoder(){}

bool GenotypeEncoder::Encode(const bcf_type& line,
		                          meta_type& meta_base,
								 block_type& block,
                           const U32* const  ppa)
{
	if(line.body->n_allele + 1 >= 32768){
		std::cerr << utility::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	meta_base.controller.biallelic        = line.body->n_allele    == 2;
	meta_base.controller.diploid          = line.gt_support.ploidy == 2;
	meta_base.controller.gt_mixed_phasing = line.gt_support.mixedPhasing;
	meta_base.controller.gt_anyMissing    = line.gt_support.hasMissing;
	meta_base.controller.gt_anyNA         = line.gt_support.hasMissing;
	meta_base.controller.gt_phase         = line.gt_support.phase;
	meta_base.controller.mixed_ploidy     = line.gt_support.hasEOV;

	// temp
	/*
	const char* const data = &line.data[line.formatID[0].l_offset];
	U32 j = 0;
	for(U32 i = 0; i < n_samples; ++i){
		const BYTE& allele1 = *reinterpret_cast<const BYTE* const>(&data[2*sizeof(BYTE)*ppa[j]]);
		const BYTE& allele2 = *reinterpret_cast<const BYTE* const>(&data[2*sizeof(BYTE)*ppa[j]+sizeof(BYTE)]);
		std::cout << (int)bcf::BCF_UNPACK_GENOTYPE(allele1) + (int)bcf::BCF_UNPACK_GENOTYPE(allele2) << '\t';
		++j;
	}
	std::cout << std::endl;
	*/

	// Assess cost and encode
	rle_helper_type cost;
	//std::cerr << meta_base.controller.biallelic << "," << line.gt_support.ploidy << "," << line.gt_support.hasEOV << std::endl;
	if(meta_base.controller.biallelic && meta_base.controller.diploid && meta_base.controller.mixed_ploidy == false){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(line, ppa);
		meta_base.controller.gt_rle = true;

		//++block.gt_rle_container;
		block.gt_support_data_container.Add((U32)cost.n_runs);
		++block.gt_support_data_container;

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(line, block.gt_rle8_container, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
			++block.gt_rle8_container;
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(line, block.gt_rle16_container, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U16;
			++block.gt_rle16_container;
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(line, block.gt_rle32_container, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U32;
			++block.gt_rle32_container;
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(line, block.gt_rle64_container, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U64;
			++block.gt_rle64_container;
			break;
		default:
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		return true;
	}
	else if(meta_base.controller.diploid) { // Case diploid n-allelic OR have EOV values
		cost = this->assessDiploidRLEnAllelic(line, ppa);

		// BCF-style cost
		U32 costBCFStyle = this->n_samples; // cost for BCF-style encoding
		if(line.body->n_allele + 1 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(line.body->n_allele + 1 < 128)   costBCFStyle *= sizeof(S16);
		else if(line.body->n_allele + 1 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width * cost.n_runs < costBCFStyle){
			meta_base.controller.gt_rle = true;
			block.gt_support_data_container.Add((U32)cost.n_runs);
			++block.gt_support_data_container;

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(line, block.gt_simple8_container, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
				++block.gt_simple8_container;
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(line, block.gt_simple16_container, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U16;
				++block.gt_simple16_container;
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(line, block.gt_simple32_container, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U32;
				++block.gt_simple32_container;
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(line, block.gt_simple64_container, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U64;
				++block.gt_simple64_container;
				break;
			default:
				std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
				return false;
			}
			return true;
		}
		// BCF style is cheaper
		else {
			//std::cerr << "BCF-style cheaper" << std::endl;
			block.gt_support_data_container.Add((U32)this->n_samples);
			++block.gt_support_data_container;

			U64 n_runs = this->n_samples;
			if(line.body->n_allele + 1 < 8){
				meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(line, block.gt_simple8_container, n_runs, ppa);
				++block.gt_simple8_container;
			}
			else if(line.body->n_allele + 1 < 128){
				meta_base.controller.gt_primtive_type = core::YON_GT_U16;
				this->EncodeDiploidBCF<U16> (line, block.gt_simple16_container, n_runs, ppa);
				++block.gt_simple16_container;
			}
			else if(line.body->n_allele + 1 < 32768){
				meta_base.controller.gt_primtive_type = core::YON_GT_U32;
				this->EncodeDiploidBCF<U32> (line, block.gt_simple64_container, n_runs, ppa);
				++block.gt_simple64_container;
			}
			else {
				std::cerr << utility::timestamp("ERROR", "ENCODER") <<
							 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
							 "Format is limited to 32768..." << std::endl;
				return false;
			}

			// Reset and recycle helper
			//this->helper.reset();
			return true;
		}
	}
	// temp
	else {
		std::cerr << "other bcf style: n_alleles: " << line.body->n_allele << ",ploidy: " << line.gt_support.ploidy << std::endl;
		block.gt_support_data_container.Add((U32)this->n_samples*line.gt_support.ploidy);
		++block.gt_support_data_container;
		++block.gt_simple8_container;

		U64 n_runs = this->n_samples*line.gt_support.ploidy;

		if(line.body->n_allele + 1 < 8)          this->EncodeBCFStyle<BYTE>(line, block.gt_simple8_container, n_runs);
		else if(line.body->n_allele + 1 < 128)   this->EncodeBCFStyle<U16> (line, block.gt_simple8_container, n_runs);
		else if(line.body->n_allele + 1 < 32768) this->EncodeBCFStyle<U32> (line, block.gt_simple8_container, n_runs);
		else {
			std::cerr << utility::timestamp("ERROR", "ENCODER") <<
						 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
						 "Format is limited to 32768..." << std::endl;
			return false;
		}
		return true;
	}
	return false;
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa) const{
	// Setup
	const BYTE ploidy = 2;
	const BYTE shift    = line.gt_support.hasMissing    ? 2 : 1;
	const BYTE add      = line.gt_support.mixedPhasing  ? 1 : 0;
	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U64 n_runs_u64  = 0; U64 run_length_u64  = 1;

	// Run limits
	const BYTE BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*shift + add)) - 1;
	const U16  U16_limit  = pow(2, 8*sizeof(U16)  - (ploidy*shift + add)) - 1;
	const U32  U32_limit  = pow(2, 8*sizeof(U32)  - (ploidy*shift + add)) - 1;
	const U64  U64_limit  = pow(2, 8*sizeof(U64)  - (ploidy*shift + add)) - 1;

	// First ref
	const char* const data = &line.data[line.formatID[0].l_offset];
	const SBYTE& allele1_2 = *reinterpret_cast<const SBYTE* const>(&data[ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2_2 = *reinterpret_cast<const SBYTE* const>(&data[ploidy*sizeof(SBYTE)*ppa[0] + sizeof(SBYTE)]);
	U32 ref = YON_PACK_GT_DIPLOID(allele2_2, allele1_2, shift, add);

	// Cycle over GT values
	U32 ppa_pos = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&data[ploidy*sizeof(SBYTE)*ppa[ppa_pos]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&data[ploidy*sizeof(SBYTE)*ppa[ppa_pos] + sizeof(SBYTE)]);
		U32 internal = YON_PACK_GT_DIPLOID(allele2, allele1, shift, add);

		// Extend or break run
		if(ref != internal){
			++n_runs_byte; run_length_byte = 0;
			++n_runs_u16;  run_length_u16  = 0;
			++n_runs_u32;  run_length_u32  = 0;
			++n_runs_u64;  run_length_u64  = 0;
			ref = internal;
		}

		// Overflow: trigger a break
		if(run_length_byte == BYTE_limit){ ++n_runs_byte; run_length_byte = 0; }
		if(run_length_u16  == U16_limit) { ++n_runs_u16;  run_length_u16  = 0; }
		if(run_length_u32  == U32_limit) { ++n_runs_u32;  run_length_u32  = 0; }
		if(run_length_u64  == U64_limit) { ++n_runs_u64;  run_length_u64  = 0; }

		// Update all counts
		++run_length_byte;
		++run_length_u16;
		++run_length_u32;
		++run_length_u64;
		++ppa_pos;
	}
	// Final runs
	++n_runs_byte;
	++n_runs_u16;
	++n_runs_u32;
	++n_runs_u64;

	// Determine best action
	U32 smallest_cost = n_runs_byte*sizeof(BYTE);
	U64 chosen_runs = n_runs_byte;
	BYTE word_width = sizeof(BYTE);
	if(n_runs_u16*sizeof(U16) < smallest_cost){ smallest_cost = n_runs_u16*sizeof(U16); word_width = sizeof(U16); chosen_runs = n_runs_u16; }
	if(n_runs_u32*sizeof(U32) < smallest_cost){ smallest_cost = n_runs_u32*sizeof(U32); word_width = sizeof(U32); chosen_runs = n_runs_u32; }
	if(n_runs_u64*sizeof(U64) < smallest_cost){ smallest_cost = n_runs_u64*sizeof(U64); word_width = sizeof(U64); chosen_runs = n_runs_u64; }

	assert(ppa_pos == n_samples);
	return(rle_helper_type(word_width, chosen_runs));
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa) const{
	const BYTE ploidy    = 2;

	// Assess RLE cost
	const BYTE shift = ceil(log2(line.body->n_allele + line.gt_support.hasMissing + line.gt_support.hasEOV + 2));
	const BYTE add   = line.gt_support.mixedPhasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size
	// Ploidy*shift_size bits for alleles and 1 bit for phase information
	S32 BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*shift + add)) - 1;
	S32  U16_limit = pow(2, 8*sizeof(U16)  - (ploidy*shift + add)) - 1;
	S64  U32_limit = pow(2, 8*sizeof(U32)  - (ploidy*shift + add)) - 1;
	U64  U64_limit = pow(2, 8*sizeof(U64)  - (ploidy*shift + add)) - 1;
	if(BYTE_limit <= 0) BYTE_limit = std::numeric_limits<S32>::max();
	if(U16_limit <= 0)  U16_limit  = std::numeric_limits<S32>::max();
	if(U32_limit <= 0)  U32_limit  = std::numeric_limits<S64>::max();

	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U64 n_runs_u64  = 0; U64 run_length_u64  = 1;

	// Setup first
	const char* const data = &line.data[line.formatID[0].l_offset];
	BYTE allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0]]);
	BYTE allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0] + sizeof(BYTE)]);
	if((allele1 >> 1) == 0){

	}
	else if(allele1 == 0x81){
		allele1 = 1;
		//std::cerr << "eov" << std::endl;
	} else {
		// Add 1 to value
		allele1 = (((allele1 >> 1) + 1) << 1) | (allele1 & 1);
	}

	if((allele2 >> 1) == 0){

	}
	else if(allele2 == 0x81){
		allele2 = 1;
		//std::cerr << "eov" << std::endl;
	} else {
		allele2 = (((allele2 >> 1) + 1) << 1) | (allele2 & 1);
	}
	U32 ref = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add);

	U32 ppa_pos = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		BYTE allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos]]);
		BYTE allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos] + sizeof(BYTE)]);
		if((allele1 >> 1) == 0){

		}
		else if(allele1 == 0x81){
			allele1 = 1;
			//std::cerr << "eov" << std::endl;
		} else {
			// Add 1 to value
			allele1 = (((allele1 >> 1) + 1) << 1) | (allele1 & 1);
		}

		if((allele2 >> 1) == 0){

		}
		else if(allele2 == 0x81){
			allele2 = 1;
			//std::cerr << "eov" << std::endl;
		} else {
			allele2 = (((allele2 >> 1) + 1) << 1) | (allele2 & 1);
		}
		const U32 internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add);

		if(ref != internal){
			ref = internal;
			++n_runs_byte; run_length_byte = 0;
			++n_runs_u16;  run_length_u16  = 0;
			++n_runs_u32;  run_length_u32  = 0;
			++n_runs_u64;  run_length_u64  = 0;
		}

		// Overflow: trigger a break
		if(run_length_byte == BYTE_limit){ ++n_runs_byte; run_length_byte = 0; }
		if(run_length_u16  == U16_limit) { ++n_runs_u16;  run_length_u16  = 0; }
		if(run_length_u32  == U32_limit) { ++n_runs_u32;  run_length_u32  = 0; }
		if(run_length_u64  == U64_limit) { ++n_runs_u64;  run_length_u64  = 0; }

		// Update all counts
		++run_length_byte;
		++run_length_u16;
		++run_length_u32;
		++run_length_u64;
		++ppa_pos;
	}
	// Final runs
	++n_runs_byte;
	++n_runs_u16;
	++n_runs_u32;
	++n_runs_u64;

	// Determine best action
	U32 smallest_cost = n_runs_byte*sizeof(BYTE);
	U64 chosen_runs = n_runs_byte;
	BYTE word_width = 1;
	if(BYTE_limit == std::numeric_limits<S32>::max()) smallest_cost = std::numeric_limits<U32>::max();
	if(n_runs_u16*sizeof(U16) < smallest_cost){ smallest_cost = n_runs_u16*sizeof(U16); word_width = sizeof(U16); chosen_runs = n_runs_u16; }
	if(n_runs_u32*sizeof(U32) < smallest_cost){ smallest_cost = n_runs_u32*sizeof(U32); word_width = sizeof(U32); chosen_runs = n_runs_u32; }
	if(n_runs_u64*sizeof(U64) < smallest_cost){ smallest_cost = n_runs_u64*sizeof(U64); word_width = sizeof(U64); chosen_runs = n_runs_u64; }

	assert(ppa_pos == n_samples);
	return(rle_helper_type(word_width, chosen_runs));
}

}
}
