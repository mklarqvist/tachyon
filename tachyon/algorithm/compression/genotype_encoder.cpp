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
							 container_type& runs,
							 container_type& simple,
                             container_type& support,
                           const U32* const  ppa)
{
	if(line.body->n_allele + 1 >= 32768){
		std::cerr << utility::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	// This is requirement for stride data
	if(support.n_entries == 0){
		support.header.controller.type              = YON_TYPE_32B;
		support.header.controller.signedness        = 0;
		support.header_stride.controller.type       = YON_TYPE_32B;
		support.header_stride.controller.signedness = 0;
	}

	meta_base.controller.biallelic        = line.body->n_allele == 2;
	meta_base.controller.diploid          = line.gt_support.ploidy == 2;
	meta_base.controller.gt_mixed_phasing = line.gt_support.mixedPhasing;
	meta_base.controller.gt_anyMissing    = line.gt_support.hasMissing;
	meta_base.controller.gt_anyNA         = line.gt_support.hasEOV;
	meta_base.controller.gt_phase         = line.gt_support.phase;
	meta_base.controller.mixed_ploidy     = line.gt_support.hasEOV;

	// Assess cost and encode
	rle_helper_type cost;
	if(line.body->n_allele == 2 && line.gt_support.ploidy == 2 && line.gt_support.hasEOV == false){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(line, ppa);
		if(!support.checkStrideSize(1))
			support.triggerMixedStride();

		support.addStride(1);
		meta_base.controller.gt_rle = true;

		++runs;
		support += (U32)cost.n_runs;
		++support;

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(line, runs, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(line, runs, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U16;
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(line, runs, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U32;
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(line, runs, ppa, cost);
			meta_base.controller.gt_primtive_type = core::YON_GT_U64;
			break;
		default:
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		return true;
	}
	else if(line.gt_support.ploidy == 2) { // Case diploid n-allelic OR have EOV values
		cost = this->assessDiploidRLEnAllelic(line, ppa);

		// BCF-style cost
		U32 costBCFStyle = this->n_samples; // cost for BCF-style encoding
		if(line.body->n_allele + 1 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(line.body->n_allele + 1 < 128)   costBCFStyle *= sizeof(S16);
		else if(line.body->n_allele + 1 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width * cost.n_runs < costBCFStyle){
			if(!support.checkStrideSize(2))
				support.triggerMixedStride();

			support.addStride(2);

			meta_base.controller.gt_rle = true;
			++simple;
			support += (U32)cost.n_runs;
			++support;

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(line, simple, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(line, simple, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U16;
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(line, simple, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U32;
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(line, simple, ppa, cost);
				meta_base.controller.gt_primtive_type = core::YON_GT_U64;
				break;
			default:
				std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
				return false;
			}
			return true;
		}
		// BCF style is cheaper
		else {
			std::cerr << "BCF-style cheaper" << std::endl;
			if(!support.checkStrideSize(3))
				support.triggerMixedStride();

			support.addStride(3);
			support += (U32)this->n_samples;
			++support;
			++simple;

			U64 n_runs = this->n_samples;
			if(line.body->n_allele + 1 < 8){
				meta_base.controller.gt_primtive_type = core::YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(line, simple, n_runs, ppa);
			}
			else if(line.body->n_allele + 1 < 128){
				meta_base.controller.gt_primtive_type = core::YON_GT_U16;
				this->EncodeDiploidBCF<U16> (line, simple, n_runs, ppa);
			}
			else if(line.body->n_allele + 1 < 32768){
				meta_base.controller.gt_primtive_type = core::YON_GT_U32;
				this->EncodeDiploidBCF<U32> (line, simple, n_runs, ppa);
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
		std::cerr << "other bcf style" << std::endl;
		if(!support.checkStrideSize(4))
			support.triggerMixedStride();

		support.addStride(4);
		support += (BYTE)line.gt_support.ploidy;
		++support;

		meta_base.controller.biallelic        = false;
		meta_base.controller.gt_mixed_phasing = true;
		meta_base.controller.gt_anyMissing    = true;
		meta_base.controller.biallelic        = line.body->n_allele == 2;
		meta_base.controller.gt_anyNA         = line.gt_support.hasEOV;
		meta_base.controller.gt_phase         = 0;
		meta_base.controller.diploid          = true;
		meta_base.controller.mixed_ploidy     = false;
		++simple;

		U64 n_runs = this->n_samples*line.gt_support.ploidy;
		if(line.body->n_allele + 1 < 8)          this->EncodeBCFStyle<BYTE>(line, simple, n_runs);
		else if(line.body->n_allele + 1 < 128)   this->EncodeBCFStyle<U16> (line, simple, n_runs);
		else if(line.body->n_allele + 1 < 32768) this->EncodeBCFStyle<U32> (line, simple, n_runs);
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
	U32 internal_buffer_offset = line.formatID[0].l_offset;
	const BYTE ploidy = 2;

	// Setup
	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U64 n_runs_u64  = 0; U64 run_length_u64  = 1;

	// First ref
	const SBYTE& allele1_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[0] + 1]);
	U32 ref = YON_PACK_GT_DIPLOID(allele2_2, allele1_2, 2, line.gt_support.mixedPhasing);

	// Run limits
	const BYTE BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*(1+line.gt_support.hasMissing)+line.gt_support.mixedPhasing)) - 1;
	const U16  U16_limit  = pow(2, 8*sizeof(U16)  - (ploidy*(1+line.gt_support.hasMissing)+line.gt_support.mixedPhasing)) - 1;
	const U32  U32_limit  = pow(2, 8*sizeof(U32)  - (ploidy*(1+line.gt_support.hasMissing)+line.gt_support.mixedPhasing)) - 1;
	const U64  U64_limit  = pow(2, 8*sizeof(U64)  - (ploidy*(1+line.gt_support.hasMissing)+line.gt_support.mixedPhasing)) - 1;

	// Cycle over GT values
	U32 j = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[j] + 1]);
		U32 internal = YON_PACK_GT_DIPLOID(allele2, allele1, 2, line.gt_support.mixedPhasing);

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

	return(rle_helper_type(word_width, chosen_runs));
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa) const{
	const BYTE ploidy    = 2;
	U32 internal_pos_rle = line.formatID[0].l_offset;

	// Assess RLE cost
	const BYTE shift     = ceil(log2(line.body->n_allele + line.gt_support.hasMissing + line.gt_support.hasEOV + 1));
	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[0] + 1]);
	U32 ref = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, line.gt_support.mixedPhasing);

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size
	// Ploidy*shift_size bits for alleles and 1 bit for phase information
	S32 BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*shift + line.gt_support.mixedPhasing + line.gt_support.hasEOV)) - 1;
	S32  U16_limit = pow(2, 8*sizeof(U16)  - (ploidy*shift + line.gt_support.mixedPhasing + line.gt_support.hasEOV)) - 1;
	S64  U32_limit = pow(2, 8*sizeof(U32)  - (ploidy*shift + line.gt_support.mixedPhasing + line.gt_support.hasEOV)) - 1;
	U64  U64_limit = pow(2, 8*sizeof(U64)  - (ploidy*shift + line.gt_support.mixedPhasing + line.gt_support.hasEOV)) - 1;
	if(BYTE_limit <= 0) BYTE_limit = std::numeric_limits<S32>::max();
	if(U16_limit <= 0)  U16_limit  = std::numeric_limits<S32>::max();
	if(U32_limit <= 0)  U32_limit  = std::numeric_limits<S64>::max();

	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U64 n_runs_u64  = 0; U64 run_length_u64  = 1;

	U32 j = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[j] + 1]);
		const U32 internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, line.gt_support.mixedPhasing);

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

	return(rle_helper_type(word_width, chosen_runs));
}

}
}
