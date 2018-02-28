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
							 const U32* const ppa){
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

	// Assess cost and encode
	rle_helper_type cost;
	if(line.body->n_allele == 2){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(line, ppa);
		if(!support.checkStrideSize(1))
			support.triggerMixedStride();

		support.addStride(1);
		meta_base.controller.gt_rle           = true;
		meta_base.controller.gt_mixed_phasing = cost.mixedPhasing;
		meta_base.controller.gt_anyMissing    = cost.hasMissing;
		meta_base.controller.biallelic        = true;
		meta_base.controller.gt_anyNA         = cost.anyNA;
		meta_base.controller.gt_phase         = cost.phased;
		meta_base.controller.diploid          = true;
		meta_base.controller.mixed_ploidy     = false;

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
	else { // Case diploid n-allelic
		cost = this->assessDiploidRLEnAllelic(line, ppa);

		meta_base.controller.biallelic        = false;
		meta_base.controller.gt_mixed_phasing = cost.mixedPhasing;
		meta_base.controller.gt_anyMissing    = cost.hasMissing;
		meta_base.controller.biallelic        = false;
		meta_base.controller.gt_anyNA         = cost.anyNA;
		meta_base.controller.gt_phase         = cost.phased;
		meta_base.controller.diploid          = true;
		meta_base.controller.mixed_ploidy     = false;

		//std::cerr << (int)cost.word_width << '\t' << cost.n_runs << '\t' << cost.word_width*cost.n_runs << '\t' << costBCFStyle << '\t' << costBCFStyle/double(cost.word_width*cost.n_runs) << std::endl;

		// BCF-style cost
		const BYTE ploidy = 2;
		U32 costBCFStyle = ploidy*this->n_samples;
		if(line.body->n_allele + 1 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(line.body->n_allele + 1 < 128)   costBCFStyle *= sizeof(S16);
		else if(line.body->n_allele + 1 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width*cost.n_runs < costBCFStyle){
			if(!support.checkStrideSize(2))
				support.triggerMixedStride();

			support.addStride(2);
			//std::cerr << line.body->POS+1 << "\t1\t0" << '\t' << (int)cost.word_width << '\t' << cost.n_runs << '\t' << (int)cost.hasMissing << '\t' << (int)cost.mixedPhasing << '\t' << cost.n_runs*cost.word_width << std::endl;

			meta_base.controller.gt_rle = true;
			++simple;
			support += (U32)cost.n_runs;
			++support;

			//std::cerr << utility::timestamp("DEBUG") << "Cost: " << cost.word_width*cost.n_runs << " @ " << cost.n_runs << '\t' << (int)cost.word_width << std::endl;

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

			// Reset and recycle helper
			//this->helper.reset();
			return true;
		}
		// BCF style is cheaper
		else {
			if(!support.checkStrideSize(3))
				support.triggerMixedStride();

			support.addStride(3);
			support += (U32)this->n_samples*2;
			++support;
			//std::cerr << line.body->POS+1 << "\t1\t1" << '\t' << (int)cost.word_width << '\t' << cost.n_runs << '\t' << (int)cost.hasMissing << '\t' << (int)cost.mixedPhasing << '\t' << costBCFStyle << std::endl;

			meta_base.controller.gt_rle = false;
			meta_base.controller.gt_primtive_type = 0;
			++simple;

			U64 n_runs = this->n_samples*2;
			if(line.body->n_allele + 1 < 8)          this->EncodeDiploidBCF<BYTE>(line, simple, n_runs, ppa);
			else if(line.body->n_allele + 1 < 128)   this->EncodeDiploidBCF<U16> (line, simple, n_runs, ppa);
			else if(line.body->n_allele + 1 < 32768) this->EncodeDiploidBCF<U32> (line, simple, n_runs, ppa);
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
	return false;
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa){
	// Assess RLE cost
	const BYTE ploidy = 2;
	U32 internal_buffer_offset = line.formatID[0].l_offset;
	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset++]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset++]);

	//\///////////////////
	// Step 1: Assessment
	//
	// Check if there are missing values
	// and/or mixed phasing
	//\///////////////////
	// Flags
	bool mixedPhase = false;
	bool anyMissing = false;
	bool anyNA      = false;

	// Set first phase for comparison
	const BYTE firstPhase = (allele2 & 1);

	// Any data is missing?
	if((allele1 >> 1) == 0 || (allele2 >> 1) == 0) anyMissing = true;

	// Cycle over GT values
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset++]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset++]);
		const BYTE& unsigned_allele1 = *reinterpret_cast<const BYTE* const>(&allele1);
		const BYTE& unsigned_allele2 = *reinterpret_cast<const BYTE* const>(&allele2);

		// Any data missing?
		if((allele1 >> 1) == 0 || (allele2 >> 1) == 0) anyMissing = true;

		// Any data NA?
		if((unsigned_allele1 >> 1) == 0x80 ||
		   (unsigned_allele2 >> 1) == 0x80 ||
		   (unsigned_allele1 >> 1) == 0x81 ||
		   (unsigned_allele2 >> 1) == 0x81) anyNA = true;

		// Any mixed phasing?
		if((allele2 & 1) != firstPhase) mixedPhase = true;
	}

	if(anyNA){
		std::cerr << "has NA" << std::endl;
		exit(1);
	}

	// Reset
	internal_buffer_offset = line.formatID[0].l_offset;

	//\///////////////////
	// Step 2
	//
	// Calculate cost / given a data type
	//\///////////////////
	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U32 n_runs_u64  = 0; U32 run_length_u64  = 1;

	// First ref
	const SBYTE& allele1_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[0] + 1]);
	U32 ref = YON_PACK_GT_DIPLOID(allele2_2, allele1_2, 2, mixedPhase);

	// Run limits
	const BYTE BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*(1+anyMissing)+mixedPhase)) - 1;
	const U16  U16_limit  = pow(2, 8*sizeof(U16)  - (ploidy*(1+anyMissing)+mixedPhase)) - 1;
	const U32  U32_limit  = pow(2, 8*sizeof(U32)  - (ploidy*(1+anyMissing)+mixedPhase)) - 1;
	const U64  U64_limit  = pow(2, 8*sizeof(U64)  - (ploidy*(1+anyMissing)+mixedPhase)) - 1;

	// Cycle over GT values
	U32 j = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_buffer_offset + ploidy*sizeof(SBYTE)*ppa[j] + 1]);
		U32 internal = YON_PACK_GT_DIPLOID(allele2, allele1, 2, mixedPhase);

		// Extend or break run
		if(ref != internal){
			//std::cerr << run_length_byte << '\t' << run_length_u16 << '\t' << run_length_u32 << '\t' << run_length_u64 << std::endl;
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

	return(rle_helper_type(word_width, chosen_runs, firstPhase, mixedPhase, anyMissing, anyNA));
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa){
	const BYTE ploidy        = 2;
	U32 internal_pos_rle     = line.formatID[0].l_offset;
	const SBYTE& pre_allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	const SBYTE& pre_allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);

	//\///////////////////
	// Step 1: Assessment
	//
	// Check if there are missing values
	// and/or mixed phasing
	//\///////////////////
	// Flags
	bool mixedPhase = false;
	bool anyMissing = false;
	bool anyNA      = false;

	// Set first phase for comparison
	const BYTE firstPhase = (pre_allele2 & 1);

	// Any data is missing?
	if((pre_allele1 >> 1) == 0 || (pre_allele2 >> 1) == 0) anyMissing = true;

	// Cycle over GT values
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		const BYTE& unsigned_allele1 = *reinterpret_cast<const BYTE* const>(&allele1);
		const BYTE& unsigned_allele2 = *reinterpret_cast<const BYTE* const>(&allele2);

		// Any data missing?
		if((allele1 >> 1) == 0 || (allele2 >> 1) == 0) anyMissing = true;

		// Any data NA?
		if((unsigned_allele1 >> 1) == 0x80 ||
		   (unsigned_allele2 >> 1) == 0x80 ||
		   (unsigned_allele1 >> 1) == 0x81 ||
		   (unsigned_allele2 >> 1) == 0x81) anyNA = true;

		// Any mixed phasing?
		if((allele2 & 1) != firstPhase) mixedPhase = true;
	}

	if(anyNA){
		std::cerr << "has NA" << std::endl;
		exit(1);
	}

	// Reset
	internal_pos_rle = line.formatID[0].l_offset;

	// Assess RLE cost
	const BYTE shift     = ceil(log2(line.body->n_allele + anyMissing + 1));
	const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[0]]);
	const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[0] + 1]);
	U32 ref = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, mixedPhase);

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size
	// Ploidy*shift_size bits for alleles and 1 bit for phase information
	S32 BYTE_limit = pow(2, 8*sizeof(BYTE) - (ploidy*shift + mixedPhase)) - 1;
	S32  U16_limit = pow(2, 8*sizeof(U16)  - (ploidy*shift + mixedPhase)) - 1;
	S64  U32_limit = pow(2, 8*sizeof(U32)  - (ploidy*shift + mixedPhase)) - 1;
	U64  U64_limit = pow(2, 8*sizeof(U64)  - (ploidy*shift + mixedPhase)) - 1;
	if(BYTE_limit <= 0) BYTE_limit = std::numeric_limits<S32>::max();
	if(U16_limit <= 0)  U16_limit  = std::numeric_limits<S32>::max();
	if(U32_limit <= 0)  U32_limit  = std::numeric_limits<S64>::max();

	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U32 n_runs_u64  = 0; U32 run_length_u64  = 1;

	U32 j = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy, ++j){
		const SBYTE& allele1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[j]]);
		const SBYTE& allele2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle + ploidy*sizeof(SBYTE)*ppa[j] + 1]);
		const U32 internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, mixedPhase);

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

	return(rle_helper_type(word_width, chosen_runs, firstPhase, mixedPhase, anyMissing, anyNA));
}

}
}
