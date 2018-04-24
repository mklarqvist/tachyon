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

bool GenotypeEncoder::Encode(const bcf_type& bcf_entry,
		                          meta_type& meta,
								 block_type& block,
                           const U32* const  ppa)
{
	if(bcf_entry.body->n_allele + 1 >= 32768){
		std::cerr << utility::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	assert(bcf_entry.body != nullptr);

	meta.controller.biallelic        = bcf_entry.body->n_allele    == 2;
	meta.controller.diploid          = bcf_entry.gt_support.ploidy == 2;
	meta.controller.gt_mixed_phasing = bcf_entry.gt_support.mixedPhasing;
	meta.controller.gt_anyMissing    = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_anyNA         = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;
	meta.controller.mixed_ploidy     = bcf_entry.gt_support.hasEOV;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;

	// Assess cost and encode
	rle_helper_type cost;
	if(meta.controller.biallelic && meta.controller.diploid && meta.controller.mixed_ploidy == false){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(bcf_entry, ppa);

		meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
		block.gt_support_data_container.Add((U32)cost.n_runs);
		++block.gt_support_data_container;

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(bcf_entry, block.gt_rle8_container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			++block.gt_rle8_container;
			++this->stats_.rle_counts[0];
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(bcf_entry, block.gt_rle16_container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U16;
			++block.gt_rle16_container;
			++this->stats_.rle_counts[1];
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(bcf_entry, block.gt_rle32_container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U32;
			++block.gt_rle32_container;
			++this->stats_.rle_counts[2];
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(bcf_entry, block.gt_rle64_container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U64;
			++block.gt_rle64_container;
			++this->stats_.rle_counts[3];
			break;
		default:
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		return true;
	}
	else if(meta.controller.diploid) { // Case diploid n-allelic OR have EOV values
		cost = this->assessDiploidRLEnAllelic(bcf_entry, ppa);

		// BCF-style cost
		U32 costBCFStyle = this->n_samples; // cost for BCF-style encoding
		if(bcf_entry.body->n_allele + 1 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(bcf_entry.body->n_allele + 1 < 128)   costBCFStyle *= sizeof(S16);
		else if(bcf_entry.body->n_allele + 1 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width * cost.n_runs < costBCFStyle){
			meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
			block.gt_support_data_container.Add((U32)cost.n_runs);
			++block.gt_support_data_container;

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(bcf_entry, block.gt_simple8_container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				++block.gt_simple8_container;
				++this->stats_.rle_simple_counts[0];
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(bcf_entry, block.gt_simple16_container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U16;
				++block.gt_simple16_container;
				++this->stats_.rle_simple_counts[1];
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(bcf_entry, block.gt_simple32_container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U32;
				++block.gt_simple32_container;
				++this->stats_.rle_simple_counts[2];
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(bcf_entry, block.gt_simple64_container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U64;
				++block.gt_simple64_container;
				++this->stats_.rle_simple_counts[3];
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

			meta.controller.gt_compression_type = YON_GT_BCF_DIPLOID;
			block.gt_support_data_container.Add((U32)this->n_samples);
			++block.gt_support_data_container;


			U64 n_runs = this->n_samples;
			if(bcf_entry.body->n_allele + 1 < 8){
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(bcf_entry, block.gt_simple8_container, n_runs, ppa);
				++block.gt_simple8_container;
				++this->stats_.diploid_bcf_counts[0];
			}
			else if(bcf_entry.body->n_allele + 1 < 128){
				meta.controller.gt_primtive_type = YON_GT_U16;
				this->EncodeDiploidBCF<U16> (bcf_entry, block.gt_simple16_container, n_runs, ppa);
				++block.gt_simple16_container;
				++this->stats_.diploid_bcf_counts[1];
			}
			else if(bcf_entry.body->n_allele + 1 < 32768){
				meta.controller.gt_primtive_type = YON_GT_U32;
				this->EncodeDiploidBCF<U32> (bcf_entry, block.gt_simple32_container, n_runs, ppa);
				++block.gt_simple32_container;
				++this->stats_.diploid_bcf_counts[2];
			}
			else {
				std::cerr << utility::timestamp("ERROR", "ENCODER") <<
							 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
							 "Format is limited to 32768..." << std::endl;
				return false;
			}
			return true;
		}
	}
	// temp
	else {
		std::cerr << "other bcf style: n_alleles: " << bcf_entry.body->n_allele << ",ploidy: " << bcf_entry.gt_support.ploidy << std::endl;
		meta.controller.gt_compression_type = YON_GT_BCF_STYLE;
		block.gt_support_data_container.Add((U32)this->n_samples*bcf_entry.gt_support.ploidy);
		++block.gt_support_data_container;

		U64 n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		if(bcf_entry.body->n_allele + 1 < 8){
			this->EncodeBCFStyle<BYTE>(bcf_entry, block.gt_simple8_container, n_runs);
			++block.gt_simple8_container;
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			++this->stats_.bcf_counts[0];
		}
		else if(bcf_entry.body->n_allele + 1 < 128){
			this->EncodeBCFStyle<U16> (bcf_entry, block.gt_simple16_container, n_runs);
			++block.gt_simple16_container;
			meta.controller.gt_primtive_type = YON_GT_U16;
			++this->stats_.bcf_counts[1];
		}
		else if(bcf_entry.body->n_allele + 1 < 32768){
			this->EncodeBCFStyle<U32> (bcf_entry, block.gt_simple32_container, n_runs);
			++block.gt_simple32_container;
			meta.controller.gt_primtive_type = YON_GT_U32;
			++this->stats_.bcf_counts[2];
		}
		else {
			std::cerr << utility::timestamp("ERROR", "ENCODER") <<
						 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
						 "Format is limited to 32768..." << std::endl;
			return false;
		}
		return true;
	}
	return false;
}

bool GenotypeEncoder::EncodeParallel(const bcf_reader_type& bcf_reader,
		                                         meta_type* meta_entries,
												block_type& block,
												 const U32* const ppa) const
{

	U32 n_threads = 28;

	GenotypeEncoderSlaveHelper* helpers = new GenotypeEncoderSlaveHelper[bcf_reader.size()];
	CalcSlave* slaves = new CalcSlave[n_threads];
	std::vector<std::thread*> threads(n_threads);
	for(U32 i = 0; i < n_threads; ++i) threads[i] = slaves[i].Start(*this, i, n_threads, bcf_reader, meta_entries, ppa, helpers);
	for(U32 i = 0; i < n_threads; ++i) threads[i]->join();
	//std::thread* t = test_slave.Start(*this, 0, 1, bcf_reader, meta_entries, ppa, helpers);
	//t->join();
	//for(U32 i = 0; i < bcf_reader.size(); ++i){
		//if(!this->EncodeParallel(bcf_reader[i], meta_entries[i], ppa, helpers[i])){
		//	std::cerr << "failed to encode" << std::endl;
		//	return false;
		//}
	//}
	// Reduce
	for(U32 i = 0; i < bcf_reader.size(); ++i) block += helpers[i];

	delete [] slaves;
	delete [] helpers;

	// Todo: outside
	//for(U32 i = thread_id; i < N; i += n_threads){
	//
	//}
	return true;
}

// Todo: pass thread with slave in it
bool GenotypeEncoder::EncodeParallel(const bcf_type& bcf_entry,
		                                  meta_type& meta,
                                   const U32* const  ppa,
						 GenotypeEncoderSlaveHelper& slave) const
{
	if(bcf_entry.body->n_allele + 1 >= 32768){
		std::cerr << utility::timestamp("ERROR", "ENCODER") <<
					 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
					 "Format is limited to 32768..." << std::endl;
		return false;
	}

	if(bcf_entry.hasGenotypes){
		meta.controller.gt_available = true;
	} else {
		meta.controller.gt_available = false;
		return true;
	}


	assert(bcf_entry.body != nullptr);

	meta.controller.biallelic        = bcf_entry.body->n_allele    == 2;
	meta.controller.diploid          = bcf_entry.gt_support.ploidy == 2;
	meta.controller.gt_mixed_phasing = bcf_entry.gt_support.mixedPhasing;
	meta.controller.gt_anyMissing    = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_anyNA         = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;
	meta.controller.mixed_ploidy     = bcf_entry.gt_support.hasEOV;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;

	U32 start_capacity = this->n_samples * 2 / 10;
	if(this->n_samples * 2 / 10 < 65536) start_capacity = 65536;
	slave.container.buffer_data_uncompressed.resize(start_capacity);

	//GenotypeEncoderSlaveHelper slave(start_capacity);

	// Assess cost and encode
	rle_helper_type cost;
	if(meta.controller.biallelic && meta.controller.diploid && meta.controller.mixed_ploidy == false){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(bcf_entry, ppa);

		slave.encoding_type = YON_GT_RLE_DIPLOID_BIALLELIC;
		slave.n_runs = cost.n_runs;

		meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
		//block.gt_support_data_container.Add((U32)cost.n_runs);
		//++block.gt_support_data_container;

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(bcf_entry, slave.container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			slave.gt_primitive = YON_GT_BYTE;
			//++block.gt_rle8_container;
			//++this->stats_.rle_counts[0];
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(bcf_entry, slave.container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U16;
			slave.gt_primitive = YON_GT_U16;
			//++block.gt_rle16_container;
			//++this->stats_.rle_counts[1];
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(bcf_entry, slave.container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U32;
			slave.gt_primitive = YON_GT_U32;
			//++block.gt_rle32_container;
			//++this->stats_.rle_counts[2];
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(bcf_entry, slave.container, ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U64;
			slave.gt_primitive = YON_GT_U64;
			//++block.gt_rle64_container;
			//++this->stats_.rle_counts[3];
			break;
		default:
			std::cerr << utility::timestamp("ERROR","ENCODER") << "Illegal word width (" << (int)cost.word_width << ")... " << std::endl;
			return false;
		}

		return true;
	}
	else if(meta.controller.diploid) { // Case diploid n-allelic OR have EOV values
		cost = this->assessDiploidRLEnAllelic(bcf_entry, ppa);

		// BCF-style cost
		U32 costBCFStyle = this->n_samples; // cost for BCF-style encoding
		if(bcf_entry.body->n_allele + 1 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(bcf_entry.body->n_allele + 1 < 128)   costBCFStyle *= sizeof(S16);
		else if(bcf_entry.body->n_allele + 1 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width * cost.n_runs < costBCFStyle){
			meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
			//block.gt_support_data_container.Add((U32)cost.n_runs);

			slave.encoding_type = YON_GT_RLE_DIPLOID_NALLELIC;
			slave.n_runs = cost.n_runs;

			//++block.gt_support_data_container;

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(bcf_entry, slave.container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				slave.gt_primitive = YON_GT_BYTE;
				//++block.gt_simple8_container;
			//	++this->stats_.rle_simple_counts[0];
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(bcf_entry, slave.container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U16;
				slave.gt_primitive = YON_GT_U16;
				//++block.gt_simple16_container;
				//++this->stats_.rle_simple_counts[1];
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(bcf_entry, slave.container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U32;
				slave.gt_primitive = YON_GT_U32;
				//++block.gt_simple32_container;
				//++this->stats_.rle_simple_counts[2];
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(bcf_entry, slave.container, ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U64;
				slave.gt_primitive = YON_GT_U64;
				//++block.gt_simple64_container;
				//++this->stats_.rle_simple_counts[3];
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

			meta.controller.gt_compression_type = YON_GT_BCF_DIPLOID;
			//block.gt_support_data_container.Add((U32)this->n_samples);

			slave.encoding_type = YON_GT_BCF_DIPLOID;
			slave.n_runs = this->n_samples;

			//++block.gt_support_data_container;


			U64 n_runs = this->n_samples;
			if(bcf_entry.body->n_allele + 1 < 8){
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				slave.gt_primitive = YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple8_container;
				//++this->stats_.diploid_bcf_counts[0];
			}
			else if(bcf_entry.body->n_allele + 1 < 128){
				meta.controller.gt_primtive_type = YON_GT_U16;
				slave.gt_primitive = YON_GT_U16;
				this->EncodeDiploidBCF<U16> (bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple16_container;
				//++this->stats_.diploid_bcf_counts[1];
			}
			else if(bcf_entry.body->n_allele + 1 < 32768){
				meta.controller.gt_primtive_type = YON_GT_U32;
				slave.gt_primitive = YON_GT_U32;
				this->EncodeDiploidBCF<U32> (bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple32_container;
				//++this->stats_.diploid_bcf_counts[2];
			}
			else {
				std::cerr << utility::timestamp("ERROR", "ENCODER") <<
							 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
							 "Format is limited to 32768..." << std::endl;
				return false;
			}
			return true;
		}
	}
	// temp
	else {
		std::cerr << "other bcf style: n_alleles: " << bcf_entry.body->n_allele << ",ploidy: " << bcf_entry.gt_support.ploidy << std::endl;
		meta.controller.gt_compression_type = YON_GT_BCF_STYLE;
		//block.gt_support_data_container.Add((U32)this->n_samples*bcf_entry.gt_support.ploidy);

		slave.encoding_type = YON_GT_BCF_STYLE;
		slave.n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		//++block.gt_support_data_container;

		U64 n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		if(bcf_entry.body->n_allele + 1 < 8){
			this->EncodeBCFStyle<BYTE>(bcf_entry, slave.container, n_runs);
			//++block.gt_simple8_container;
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			slave.gt_primitive = YON_GT_BYTE;
			//++this->stats_.bcf_counts[0];
		}
		else if(bcf_entry.body->n_allele + 1 < 128){
			this->EncodeBCFStyle<U16> (bcf_entry, slave.container, n_runs);
			//++block.gt_simple16_container;
			meta.controller.gt_primtive_type = YON_GT_U16;
			slave.gt_primitive = YON_GT_U16;
			//++this->stats_.bcf_counts[1];
		}
		else if(bcf_entry.body->n_allele + 1 < 32768){
			this->EncodeBCFStyle<U32> (bcf_entry, slave.container, n_runs);
			//++block.gt_simple32_container;
			meta.controller.gt_primtive_type = YON_GT_U32;
			slave.gt_primitive = YON_GT_U32;
			//++this->stats_.bcf_counts[2];
		}
		else {
			std::cerr << utility::timestamp("ERROR", "ENCODER") <<
						 "Illegal number of alleles (" << bcf_entry.body->n_allele + 1 << "). "
						 "Format is limited to 32768..." << std::endl;
			return false;
		}
		return true;
	}
	return false;
}

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEBiallelic(const bcf_type& bcf_entry, const U32* const ppa) const{
	// Setup
	const BYTE ploidy = 2;
	const BYTE shift    = bcf_entry.gt_support.hasMissing    ? 2 : 1;
	const BYTE add      = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;
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
	const char* const data = &bcf_entry.data[bcf_entry.formatID[0].l_offset];
	const BYTE& allele1_2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0]]);
	const BYTE& allele2_2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0] + sizeof(BYTE)]);
	U32 ref = YON_PACK_GT_DIPLOID(allele2_2, allele1_2, shift, add);

	// Cycle over GT values
	U32 ppa_pos = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		const BYTE& allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos]]);
		const BYTE& allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos] + sizeof(BYTE)]);
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

const GenotypeEncoder::rle_helper_type GenotypeEncoder::assessDiploidRLEnAllelic(const bcf_type& bcf_entry, const U32* const ppa) const{
	const BYTE ploidy    = 2;

	// Assess RLE cost
	const BYTE shift = ceil(log2(bcf_entry.body->n_allele + bcf_entry.gt_support.hasMissing + bcf_entry.gt_support.hasEOV + 1));
	const BYTE add   = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;

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
	const char* const data = &bcf_entry.data[bcf_entry.formatID[0].l_offset];
	BYTE allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0]]);
	BYTE allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0] + sizeof(BYTE)]);
	if((allele1 >> 1) == 0){
		allele1 = 0;
	}
	else if(allele1 == 0x81){
		allele1 = 1;
		//std::cerr << "eov" << std::endl;
	} else {
		// Add 1 to value
		allele1 = (((allele1 >> 1) + bcf_entry.gt_support.hasEOV) << 1) | (allele1 & 1);
	}

	if((allele2 >> 1) == 0){
		allele2 = 0;
	}
	else if(allele2 == 0x81){
		allele2 = 1;
		//std::cerr << "eov" << std::endl;
	} else {
		allele2 = (((allele2 >> 1) + bcf_entry.gt_support.hasEOV) << 1) | (allele2 & 1);
	}
	U32 ref = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add);

	U32 ppa_pos = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		BYTE allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos]]);
		BYTE allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos] + sizeof(BYTE)]);
		if((allele1 >> 1) == 0){
			allele1 = 0;
		}
		else if(allele1 == 0x81){
			allele1 = 1;
			//std::cerr << "eov" << std::endl;
		} else {
			// Add 1 to value
			allele1 = (((allele1 >> 1) + bcf_entry.gt_support.hasEOV) << 1) | (allele1 & 1);
		}

		if((allele2 >> 1) == 0){
			allele2 = 0;
		}
		else if(allele2 == 0x81){
			allele2 = 1;
			//std::cerr << "eov" << std::endl;
		} else {
			allele2 = (((allele2 >> 1) + bcf_entry.gt_support.hasEOV) << 1) | (allele2 & 1);
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
