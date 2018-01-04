#include "EncoderGenotypesRLE.h"

namespace Tachyon{
namespace Encoding{

EncoderGenotypesRLE::EncoderGenotypesRLE() :
	n_samples(0)
{
}

EncoderGenotypesRLE::EncoderGenotypesRLE(const U64 samples) :
	n_samples(samples)
	//helper(samples)
{
}

EncoderGenotypesRLE::~EncoderGenotypesRLE(){}

bool EncoderGenotypesRLE::Encode(const bcf_type& line, meta_type& meta_base, container_type& runs, container_type& simple, container_type& support, U64& n_runs, const U32* const ppa){
	if(line.body->n_allele + 1 >= 32768){
		std::cerr << Helpers::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << line.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	// Assess cost and encode
	rle_helper_type cost;
	if(line.body->n_allele == 2){
		cost = this->assessDiploidRLEBiallelic(line, ppa);
		//meta_base.virtual_offset_gt        = runs.buffer_data_uncompressed.pointer; // absolute position at start of stream
		support += (U32)0;
		meta_base.controller.rle           = true;
		meta_base.controller.mixed_phasing = cost.mixedPhasing;
		meta_base.controller.anyMissing    = cost.hasMissing;
		meta_base.controller.biallelic     = true;
		++runs.n_entries;

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(line, runs, n_runs, ppa, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 0;
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(line, runs, n_runs, ppa, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 1;
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(line, runs, n_runs, ppa, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 2;
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(line, runs, n_runs, ppa, cost.hasMissing, cost.mixedPhasing);
			meta_base.controller.rle_type = 3;
			break;
		default:
			std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		//meta_base.AF = float(this->helper.countsAlleles[1]) / (this->helper.countsAlleles[0] + this->helper.countsAlleles[1]);

		// Reset and recycle helper
		//this->helper.reset();

		return true;
	}
	else {
		cost = this->assessDiploidRLEnAllelic(line, ppa);
		U32 costBCFStyle = this->n_samples;
		if(line.body->n_allele + 1 < 8)          costBCFStyle *= 1;
		else if(line.body->n_allele + 1 < 128)   costBCFStyle *= 2;
		else if(line.body->n_allele + 1 < 32768) costBCFStyle *= 4;

		support += (U32)1;
		//meta_base.virtual_offset_gt        = simple.buffer_data_uncompressed.pointer; // absolute position at start of stream
		meta_base.controller.biallelic     = false;
		//meta_base.AF                       = 0; // AF needs to be looked up in cold store
		meta_base.controller.mixed_phasing = cost.mixedPhasing;
		meta_base.controller.anyMissing    = cost.hasMissing;
		meta_base.controller.simple        = 1;

		//std::cerr << (int)cost.word_width << '\t' << cost.n_runs << '\t' << cost.word_width*cost.n_runs << '\t' << costBCFStyle << '\t' << costBCFStyle/double(cost.word_width*cost.n_runs) << std::endl;

		// RLE is cheaper
		if(cost.word_width*cost.n_runs < costBCFStyle){
			//std::cerr << line.body->POS+1 << "\t1\t0" << '\t' << (int)cost.word_width << '\t' << cost.n_runs << '\t' << (int)cost.hasMissing << '\t' << (int)cost.mixedPhasing << '\t' << cost.n_runs*cost.word_width << std::endl;

			meta_base.controller.rle = true;
			++simple.n_entries;

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(line, simple, n_runs, ppa);
				meta_base.controller.rle_type = 0;
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(line, simple, n_runs, ppa);
				meta_base.controller.rle_type = 1;
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(line, simple, n_runs, ppa);
				meta_base.controller.rle_type = 2;
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(line, simple, n_runs, ppa);
				meta_base.controller.rle_type = 3;
				break;
			default:
				std::cerr << Helpers::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
				return false;
			}

			// Reset and recycle helper
			//this->helper.reset();
			return true;
		}
		// BCF style is cheaper
		else {
			//std::cerr << line.body->POS+1 << "\t1\t1" << '\t' << (int)cost.word_width << '\t' << cost.n_runs << '\t' << (int)cost.hasMissing << '\t' << (int)cost.mixedPhasing << '\t' << costBCFStyle << std::endl;

			meta_base.controller.rle      = false;
			meta_base.controller.rle_type = 0;
			++simple.n_entries;

			if(line.body->n_allele + 1 < 8)          this->EncodeDiploidSimple<BYTE>(line, simple, n_runs, ppa);
			else if(line.body->n_allele + 1 < 128)   this->EncodeDiploidSimple<U16> (line, simple, n_runs, ppa);
			else if(line.body->n_allele + 1 < 32768) this->EncodeDiploidSimple<U32> (line, simple, n_runs, ppa);
			else {
				std::cerr << Helpers::timestamp("ERROR", "ENCODER") <<
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

const EncoderGenotypesRLE::rle_helper_type EncoderGenotypesRLE::assessDiploidRLEBiallelic(const bcf_type& line, const U32* const ppa){
	// Assess RLE cost
	U32 internal_pos_rle = line.p_genotypes;
	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);

	/////////////////////
	// Step 1: Assessment
	//
	// Check if there are missing values
	// and/or mixed phasing
	/////////////////////
	// Flags
	bool mixedPhase = false;
	bool anyMissing = false;

	// Set first phase for comparison
	const BYTE firstPhase = (fmt_type_value2 & 1);

	// Any data is missing?
	if((fmt_type_value1 >> 1) == 0 || (fmt_type_value2 >> 1) == 0) anyMissing = true;

	// Cycle over GT values
	for(U32 i = 2; i < this->n_samples * 2; i += 2){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle++]);

		// Any data is missing?
		if((fmt_type_value1 >> 1) == 0 || (fmt_type_value2 >> 1) == 0) anyMissing = true;
		// Is there mixed phasing?
		if((fmt_type_value2 & 1) != firstPhase) mixedPhase = true;
	}

	// Reset
	internal_pos_rle = line.p_genotypes;

	/////////////////////
	// Step 2
	//
	// Calculate cost / given a data type
	/////////////////////
	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U32 n_runs_u64  = 0; U32 run_length_u64  = 1;

	// First ref
	const SBYTE& fmt_type_value1_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[0]]);
	const SBYTE& fmt_type_value2_2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[0]+1]);
	U32 ref = PACK_RLE_BIALLELIC(fmt_type_value2_2, fmt_type_value1_2, 2, 1);

	// Run limits
	const BYTE BYTE_limit = pow(2, 8*sizeof(BYTE) - (2*(1+anyMissing)+mixedPhase)) - 1;
	const U16  U16_limit  = pow(2, 8*sizeof(U16)  - (2*(1+anyMissing)+mixedPhase)) - 1;
	const U32  U32_limit  = pow(2, 8*sizeof(U32)  - (2*(1+anyMissing)+mixedPhase)) - 1;
	const U64  U64_limit  = pow(2, 8*sizeof(U64)  - (2*(1+anyMissing)+mixedPhase)) - 1;

	// Cycle over GT values
	U32 j = 1;
	for(U32 i = 2; i < this->n_samples * 2; i += 2, ++j){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[j]]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[j]+1]);
		U32 internal = PACK_RLE_BIALLELIC(fmt_type_value2, fmt_type_value1, 2, 1);

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
	BYTE word_width = 1;
	if(n_runs_u16*sizeof(U16) < smallest_cost){ smallest_cost = n_runs_u16*sizeof(U16); word_width = 2; chosen_runs = n_runs_u16; }
	if(n_runs_u32*sizeof(U32) < smallest_cost){ smallest_cost = n_runs_u32*sizeof(U32); word_width = 4; chosen_runs = n_runs_u32; }
	if(n_runs_u64*sizeof(U64) < smallest_cost){ smallest_cost = n_runs_u64*sizeof(U64); word_width = 8; chosen_runs = n_runs_u64; }

	return(rle_helper_type(word_width, chosen_runs, mixedPhase, anyMissing));
}

const EncoderGenotypesRLE::rle_helper_type EncoderGenotypesRLE::assessDiploidRLEnAllelic(const bcf_type& line, const U32* const ppa){
	// Assess RLE cost
	const BYTE shift_size = ceil(log2(line.body->n_allele + 1));
	U32 internal_pos_rle = line.p_genotypes;
	const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[0]]);
	const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[0]+1]);
	U32 ref = PACK_RLE_SIMPLE(fmt_type_value2, fmt_type_value1, shift_size);

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size
	S32 BYTE_limit  = pow(2, 8*sizeof(BYTE) - (2*(shift_size)+1)) - 1;
	S32  U16_limit  = pow(2, 8*sizeof(U16)  - (2*(shift_size)+1)) - 1;
	S32  U32_limit  = pow(2, 8*sizeof(U32)  - (2*(shift_size)+1)) - 1;
	U64  U64_limit  = pow(2, 8*sizeof(U64)  - (2*(shift_size)+1)) - 1;
	if(BYTE_limit < 0) BYTE_limit = std::numeric_limits<S32>::max();
	if(U16_limit < 0)  U16_limit  = std::numeric_limits<S32>::max();
	if(U32_limit < 0)  U32_limit  = std::numeric_limits<S32>::max();

	U32 n_runs_byte = 0; U32 run_length_byte = 1;
	U32 n_runs_u16  = 0; U32 run_length_u16  = 1;
	U32 n_runs_u32  = 0; U32 run_length_u32  = 1;
	U32 n_runs_u64  = 0; U32 run_length_u64  = 1;

	U32 j = 1;
	for(U32 i = 2; i < this->n_samples * 2; i += 2, ++j){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[j]]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&line.data[internal_pos_rle+2*ppa[j]+1]);
		U32 internal = PACK_RLE_SIMPLE(fmt_type_value2, fmt_type_value1, shift_size);

		if(ref != internal){
			ref = internal;
			++n_runs_byte; run_length_byte = 0;
			++n_runs_u16;  run_length_u16  = 0;
			++n_runs_u32;  run_length_u32  = 0;
			++n_runs_u64;  run_length_u64  = 0;
		}

		// Overflow: trigger a break
		if(run_length_byte == BYTE_limit){ ++n_runs_byte; run_length_byte = 0; }
		if(run_length_u16  == U16_limit) { ++n_runs_u16; run_length_u16   = 0; }
		if(run_length_u32  == U32_limit) { ++n_runs_u32; run_length_u32   = 0; }
		if(run_length_u64  == U64_limit) { ++n_runs_u64; run_length_u64   = 0; }

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
	if(n_runs_u16*sizeof(U16) < smallest_cost && U16_limit != std::numeric_limits<S32>::max()){ smallest_cost = n_runs_u16*sizeof(U16); word_width = 2; chosen_runs = n_runs_u16; }
	if(n_runs_u32*sizeof(U32) < smallest_cost && U32_limit != std::numeric_limits<S32>::max()){ smallest_cost = n_runs_u32*sizeof(U32); word_width = 4; chosen_runs = n_runs_u32; }
	if(n_runs_u64*sizeof(U64) < smallest_cost){ smallest_cost = n_runs_u64*sizeof(U64); word_width = 8; chosen_runs = n_runs_u64; }

	//std::cerr << (int)word_width << '\t' << chosen_runs << '\t' << word_width*chosen_runs << "\t\t" << n_runs_byte*sizeof(BYTE) << '\t' << n_runs_u16*sizeof(U16) << '\t' << n_runs_u32*sizeof(U32) << '\t' << n_runs_u64*sizeof(U64) << std::endl;
	//std::cerr << BYTE_limit << '\t' << U16_limit << '\t'  << U32_limit << '\t' << U64_limit << std::endl;

	return(rle_helper_type(word_width, chosen_runs, true, true));
}

}
}
