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

bool GenotypeEncoder::AssessDiploidBiallelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
	assert(entry->d.fmt[0].n == 2);
	assert(entry->n_allele == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_samples * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	// Track all possible outcomes.
	// 1: BYTE + Permuted
	// 2: U16  + Permuted
	// 3: U32  + Permuted
	// 4: U64  + Permuted
	// 5: BYTE + No permutation
	// 6: U16  + No permutation
	// 7: U32  + No permutation
	// 8: U64  + No permutation
	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	memset(l_runs, 1, sizeof(uint64_t)*8);
	memset(n_runs, 0, sizeof(uint64_t)*8);

	// 1 + hasMissing + hasMixedPhasing
	const BYTE shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const BYTE add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	uint64_t limits[4];
	limits[0] = pow(2, 8*sizeof(BYTE) - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(U16)  - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(U32)  - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(U64)  - (base_ploidy*shift + add)) - 1;

	U32 rle_current_ref     = YON_PACK_GT_DIPLOID(gt[0], gt[1], shift, add);
	U32 rle_ppa_current_ref = YON_PACK_GT_DIPLOID(gt[permutation_array[0] * sizeof(int8_t) * base_ploidy],
	                                              gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1],
												  shift, add);

	// Keep track of the linear offset in the genotype
	// data stream. The permuted offset is computed directly
	// and does not need to be tracked.
	uint32_t l_gt_offset = 2;

	// Iterate over all available samples.
	for(U32 i = 1; i < this->n_samples; ++i, l_gt_offset += 2){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		U32 rle_current     = YON_PACK_GT_DIPLOID(gt[l_gt_offset], gt[l_gt_offset + 1], shift, add);
		U32 rle_ppa_current = YON_PACK_GT_DIPLOID(gt_ppa_target[0], gt_ppa_target[1], shift, add);

		if(rle_current != rle_current_ref){
			++n_runs[4]; ++n_runs[5]; ++n_runs[6]; ++n_runs[7];
			memset(&l_runs[4], 0, sizeof(uint64_t)*4);
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for(U32 k = 4; k < 8; ++k){
			if(l_runs[k] == limits[k-4]){
				++n_runs[k]; l_runs[k] = 0;
			}
		}

		++l_runs[4]; ++l_runs[5]; ++l_runs[6]; ++l_runs[7];


		if(rle_ppa_current != rle_ppa_current_ref){
			++n_runs[0]; ++n_runs[1]; ++n_runs[2]; ++n_runs[3];
			memset(&l_runs[0], 0, sizeof(uint64_t)*4);
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for(U32 k = 0; k < 4; ++k){
			if(l_runs[k] == limits[k]){
				++n_runs[k]; l_runs[k] = 0;
			}
		}
		++l_runs[0]; ++l_runs[1]; ++l_runs[2]; ++l_runs[3];

	}
	assert(l_gt_offset == l_gt);
	for(U32 k = 0; k < 8; ++k) ++n_runs[k];

	std::cout << entry->pos + 1 << "\tR";
	for(U32 i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for(U32 i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;

	return true;
}

bool GenotypeEncoder::AssessDiploidMultiAllelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
	assert(entry->d.fmt[0].n == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_samples * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	// Track all possible outcomes.
	// 1: BYTE + Permuted
	// 2: U16  + Permuted
	// 3: U32  + Permuted
	// 4: U64  + Permuted
	// 5: BYTE + No permutation
	// 6: U16  + No permutation
	// 7: U32  + No permutation
	// 8: U64  + No permutation
	int64_t n_runs[8]; // Number of runs.
	int64_t l_runs[8]; // Current run length.
	memset(l_runs, 1, sizeof(int64_t)*8);
	memset(n_runs, 0, sizeof(int64_t)*8);

	// Assess RLE cost
	const BYTE shift = ceil(log2(entry->n_allele + 2 + 1));
	const BYTE add   = gt_summary.n_vector_end  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size.
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
	int64_t limits[4];
	limits[0] = pow(2, 8*sizeof(BYTE) - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(U16)  - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(U32)  - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(U64)  - (base_ploidy*shift + add)) - 1;
	if(limits[0] <= 0) limits[0] = std::numeric_limits<int64_t>::max();
	if(limits[1] <= 0) limits[1] = std::numeric_limits<int64_t>::max();
	if(limits[2] <= 0) limits[2] = std::numeric_limits<int64_t>::max();

	uint8_t gt_remap[256];
	memset(gt_remap, 256, 255);
	for(U32 i = 0; i <= entry->n_allele; ++i){
		gt_remap[i << 1]       = (i << 1) + 1;
		gt_remap[(i << 1) + 1] = (i << 1) + 2;
	}
	gt_remap[0] = 0;
	gt_remap[129] = 1;

	U32 rle_current_ref     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[0]], gt_remap[gt[1]], shift, add, gt_remap[gt[1]]);
	U32 rle_ppa_current_ref = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy]],
												           gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]],
												           shift, add,
												           gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]]);

	// Keep track of the linear offset in the genotype
	// data stream. The permuted offset is computed directly
	// and does not need to be tracked.
	uint32_t l_gt_offset = base_ploidy;

	//std::cerr << entry->pos + 1 << "\t" << (gt[0]>>1) << "|" << (gt[1]>>1);

	// Iterate over all available samples.
	for(U32 i = 1; i < this->n_samples; ++i, l_gt_offset += base_ploidy){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		//std::cerr << "\t" << (gt[l_gt_offset] >> 1) << "|" << (gt[l_gt_offset+1] >> 1);
		U32 rle_current     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[l_gt_offset]], gt_remap[gt[l_gt_offset+1]], shift, add, gt_remap[gt[1]]);
		U32 rle_ppa_current = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt_ppa_target[0]],
		                                                   gt_remap[gt_ppa_target[1]],
												           shift, add,
														   gt_remap[gt_ppa_target[1]]);

		assert(gt[l_gt_offset] != 255);
		assert(gt[l_gt_offset+1] != 255);
		assert(gt_remap[gt_ppa_target[0]] != 255);
		assert(gt_remap[gt_ppa_target[1]] != 255);

		/*
		std::cerr << (U32)gt[l_gt_offset] << "->" << (U32)gt_remap[gt[l_gt_offset]] << std::endl;
		std::cerr << (U32)gt[l_gt_offset+1] << "->" << (U32)gt_remap[gt[l_gt_offset+1]] << std::endl;
		if(gt[l_gt_offset] == 129){
			std::cerr << (int)gt[l_gt_offset] << "==129 -> " << (U32)gt_remap[gt[l_gt_offset]] << std::endl;
		}
		if(gt[l_gt_offset+1] == 129){
			std::cerr << (int)gt[l_gt_offset+1] << "==129 -> " << (U32)gt_remap[gt[l_gt_offset+1]] << std::endl;
		}
		*/


		if(rle_current != rle_current_ref){
			++n_runs[4]; ++n_runs[5]; ++n_runs[6]; ++n_runs[7];
			memset(&l_runs[4], 0, sizeof(uint64_t)*4);
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for(U32 k = 4; k < 8; ++k){
			if(l_runs[k] == limits[k-4]){
				++n_runs[k]; l_runs[k] = 0;
			}
		}

		++l_runs[4]; ++l_runs[5]; ++l_runs[6]; ++l_runs[7];


		if(rle_ppa_current != rle_ppa_current_ref){
			++n_runs[0]; ++n_runs[1]; ++n_runs[2]; ++n_runs[3];
			memset(&l_runs[0], 0, sizeof(uint64_t)*4);
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for(U32 k = 0; k < 4; ++k){
			if(l_runs[k] == limits[k]){
				++n_runs[k]; l_runs[k] = 0;
			}
		}
		++l_runs[0]; ++l_runs[1]; ++l_runs[2]; ++l_runs[3];

	}
	//std::cerr << std::endl;
	assert(l_gt_offset == l_gt);
	for(U32 k = 0; k < 8; ++k) ++n_runs[k];

	std::cout << entry->pos + 1 << "\tM";
	for(U32 i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for(U32 i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;

	return true;
}

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

	if(bcf_entry.hasGenotypes){
		meta.controller.gt_available = true;
	} else {
		meta.controller.gt_available = false;
		return true;
	}

	// Assess cost and encode
	rle_helper_type cost;
	if(meta.controller.biallelic && meta.controller.diploid && meta.controller.mixed_ploidy == false){ // Case diploid and biallelic
		cost = this->assessDiploidRLEBiallelic(bcf_entry, ppa);

		meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
		block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)cost.n_runs);
		++block.base_containers[YON_BLK_GT_SUPPORT];

		switch(cost.word_width){
		case 1:
			this->EncodeDiploidRLEBiallelic<BYTE>(bcf_entry, block.base_containers[YON_BLK_GT_INT8], ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			++block.base_containers[YON_BLK_GT_INT8];
			++this->stats_.rle_counts[0];
			break;
		case 2:
			this->EncodeDiploidRLEBiallelic<U16>(bcf_entry, block.base_containers[YON_BLK_GT_INT16], ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U16;
			++block.base_containers[YON_BLK_GT_INT8];
			++this->stats_.rle_counts[1];
			break;
		case 4:
			this->EncodeDiploidRLEBiallelic<U32>(bcf_entry, block.base_containers[YON_BLK_GT_INT32], ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U32;
			++block.base_containers[YON_BLK_GT_INT32];
			++this->stats_.rle_counts[2];
			break;
		case 8:
			this->EncodeDiploidRLEBiallelic<U64>(bcf_entry, block.base_containers[YON_BLK_GT_INT64], ppa, cost);
			meta.controller.gt_primtive_type = YON_GT_U64;
			++block.base_containers[YON_BLK_GT_INT64];
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
		if(bcf_entry.body->n_allele + 2 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(bcf_entry.body->n_allele + 2 < 128)   costBCFStyle *= sizeof(S16);
		else if(bcf_entry.body->n_allele + 2 < 32768) costBCFStyle *= sizeof(S32);

		// RLE is cheaper
		if(cost.word_width * cost.n_runs < costBCFStyle){
			meta.controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
			block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)cost.n_runs);
			++block.base_containers[YON_BLK_GT_SUPPORT];

			switch(cost.word_width){
			case 1:
				this->EncodeDiploidRLEnAllelic<BYTE>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT8], ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				++block.base_containers[YON_BLK_GT_S_INT8];
				++this->stats_.rle_simple_counts[0];
				break;
			case 2:
				this->EncodeDiploidRLEnAllelic<U16>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT16], ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U16;
				++block.base_containers[YON_BLK_GT_S_INT16];
				++this->stats_.rle_simple_counts[1];
				break;
			case 4:
				this->EncodeDiploidRLEnAllelic<U32>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT32], ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U32;
				++block.base_containers[YON_BLK_GT_S_INT32];
				++this->stats_.rle_simple_counts[2];
				break;
			case 8:
				this->EncodeDiploidRLEnAllelic<U64>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT64], ppa, cost);
				meta.controller.gt_primtive_type = YON_GT_U64;
				++block.base_containers[YON_BLK_GT_S_INT64];
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
			block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)this->n_samples);
			++block.base_containers[YON_BLK_GT_SUPPORT];


			U64 n_runs = this->n_samples;
			if(bcf_entry.body->n_allele + 2 < 8){
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT8], n_runs, ppa);
				++block.base_containers[YON_BLK_GT_S_INT8];
				++this->stats_.diploid_bcf_counts[0];
			}
			else if(bcf_entry.body->n_allele + 2 < 128){
				meta.controller.gt_primtive_type = YON_GT_U16;
				this->EncodeDiploidBCF<U16> (bcf_entry, block.base_containers[YON_BLK_GT_S_INT16], n_runs, ppa);
				++block.base_containers[YON_BLK_GT_S_INT16];
				++this->stats_.diploid_bcf_counts[1];
			}
			else if(bcf_entry.body->n_allele + 2 < 32768){
				meta.controller.gt_primtive_type = YON_GT_U32;
				this->EncodeDiploidBCF<U32> (bcf_entry, block.base_containers[YON_BLK_GT_S_INT32], n_runs, ppa);
				++block.base_containers[YON_BLK_GT_S_INT32];
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
		block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)this->n_samples*bcf_entry.gt_support.ploidy);
		++block.base_containers[YON_BLK_GT_SUPPORT];

		U64 n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		if(bcf_entry.body->n_allele + 2 < 8){
			this->EncodeBCFStyle<BYTE>(bcf_entry, block.base_containers[YON_BLK_GT_S_INT8], n_runs);
			++block.base_containers[YON_BLK_GT_S_INT8];
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			++this->stats_.bcf_counts[0];
		}
		else if(bcf_entry.body->n_allele + 2 < 128){
			this->EncodeBCFStyle<U16> (bcf_entry, block.base_containers[YON_BLK_GT_S_INT16], n_runs);
			++block.base_containers[YON_BLK_GT_S_INT16];
			meta.controller.gt_primtive_type = YON_GT_U16;
			++this->stats_.bcf_counts[1];
		}
		else if(bcf_entry.body->n_allele + 2 < 32768){
			this->EncodeBCFStyle<U32> (bcf_entry, block.base_containers[YON_BLK_GT_S_INT32], n_runs);
			++block.base_containers[YON_BLK_GT_S_INT32];
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
												 const U32* const ppa,
												 const U32 n_threads)
{
	GenotypeEncoderSlaveHelper* helpers = new GenotypeEncoderSlaveHelper[bcf_reader.size()];
	CalcSlave* slaves = new CalcSlave[n_threads];
	std::vector<std::thread*> threads(n_threads);
	for(U32 i = 0; i < n_threads; ++i) threads[i] = slaves[i].Start(*this, i, n_threads, bcf_reader, meta_entries, ppa, helpers);
	for(U32 i = 0; i < n_threads; ++i) threads[i]->join();
	for(U32 i = 0; i < bcf_reader.size(); ++i) {
		block += helpers[i];
		this->updateStatistics(helpers[i]);
	}

	delete [] slaves;
	delete [] helpers;

	return true;
}

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

	assert(bcf_entry.body != nullptr);

	meta.controller.biallelic        = bcf_entry.body->n_allele    == 2;
	meta.controller.diploid          = bcf_entry.gt_support.ploidy == 2;
	meta.controller.gt_mixed_phasing = bcf_entry.gt_support.mixedPhasing;
	meta.controller.gt_anyMissing    = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_anyNA         = bcf_entry.gt_support.hasMissing;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;
	meta.controller.mixed_ploidy     = bcf_entry.gt_support.hasEOV;
	meta.controller.gt_phase         = bcf_entry.gt_support.phase;

	if(bcf_entry.hasGenotypes){
		meta.controller.gt_available = true;
	} else {
		meta.controller.gt_available = false;
		return true;
	}

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
		if(bcf_entry.body->n_allele + 2 < 8)          costBCFStyle *= sizeof(SBYTE);
		else if(bcf_entry.body->n_allele + 2 < 128)   costBCFStyle *= sizeof(S16);
		else if(bcf_entry.body->n_allele + 2 < 32768) costBCFStyle *= sizeof(S32);

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
			if(bcf_entry.body->n_allele + 2 < 8){
				meta.controller.gt_primtive_type = YON_GT_BYTE;
				slave.gt_primitive = YON_GT_BYTE;
				this->EncodeDiploidBCF<BYTE>(bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple8_container;
				//++this->stats_.diploid_bcf_counts[0];
			}
			else if(bcf_entry.body->n_allele + 2 < 128){
				meta.controller.gt_primtive_type = YON_GT_U16;
				slave.gt_primitive = YON_GT_U16;
				this->EncodeDiploidBCF<U16>(bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple16_container;
				//++this->stats_.diploid_bcf_counts[1];
			}
			else if(bcf_entry.body->n_allele + 2 < 32768){
				meta.controller.gt_primtive_type = YON_GT_U32;
				slave.gt_primitive = YON_GT_U32;
				this->EncodeDiploidBCF<U32>(bcf_entry, slave.container, n_runs, ppa);
				//++block.gt_simple32_container;
				//++this->stats_.diploid_bcf_counts[2];
			}
			else {
				std::cerr << utility::timestamp("ERROR", "ENCODER") <<
							 "Illegal number of alleles (" << bcf_entry.body->n_allele + 2 << "). "
							 "Format is limited to 32768..." << std::endl;
				return false;
			}
			return true;
		}
	}
	else {
		//std::cerr << "other bcf style: n_alleles: " << bcf_entry.body->n_allele << ",ploidy: " << bcf_entry.gt_support.ploidy << std::endl;
		meta.controller.gt_compression_type = YON_GT_BCF_STYLE;
		//block.gt_support_data_container.Add((U32)this->n_samples*bcf_entry.gt_support.ploidy);

		slave.encoding_type = YON_GT_BCF_STYLE;
		slave.n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		//++block.gt_support_data_container;

		U64 n_runs = this->n_samples*bcf_entry.gt_support.ploidy;

		if(bcf_entry.body->n_allele + 2 < 8){
			this->EncodeBCFStyle<BYTE>(bcf_entry, slave.container, n_runs);
			//++block.gt_simple8_container;
			meta.controller.gt_primtive_type = YON_GT_BYTE;
			slave.gt_primitive = YON_GT_BYTE;
			//++this->stats_.bcf_counts[0];
		}
		else if(bcf_entry.body->n_allele + 2 < 128){
			this->EncodeBCFStyle<U16> (bcf_entry, slave.container, n_runs);
			//++block.gt_simple16_container;
			meta.controller.gt_primtive_type = YON_GT_U16;
			slave.gt_primitive = YON_GT_U16;
			//++this->stats_.bcf_counts[1];
		}
		else if(bcf_entry.body->n_allele + 2 < 32768){
			this->EncodeBCFStyle<U32> (bcf_entry, slave.container, n_runs);
			//++block.gt_simple32_container;
			meta.controller.gt_primtive_type = YON_GT_U32;
			slave.gt_primitive = YON_GT_U32;
			//++this->stats_.bcf_counts[2];
		}
		else {
			std::cerr << utility::timestamp("ERROR", "ENCODER") <<
						 "Illegal number of alleles (" << bcf_entry.body->n_allele + 2 << "). "
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
	// 1 + hasMissing + hasMixedPhasing;
	const BYTE shift  = bcf_entry.gt_support.hasMissing    ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const BYTE add    = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;
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
	const BYTE shift = ceil(log2(bcf_entry.body->n_allele + 2 + 1));
	const BYTE add   = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
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
	const bool phase = allele2 & 1;
	if((allele1 >> 1) == 0)  allele1 = 0;
	else if(allele1 == 0x81) allele1 = 1;
	else allele1 = (allele1 >> 1) + 1;

	if((allele2 >> 1) == 0)  allele2 = 0;
	else if(allele2 == 0x81) allele2 = 1;
	else allele2 = (allele2 >> 1) + 1;
	U32 ref = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add, phase);

	U32 ppa_pos = 1;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		BYTE allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos]]);
		BYTE allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos] + sizeof(BYTE)]);
		const bool phase = allele2 & 1;

		if((allele1 >> 1) == 0)  allele1 = 0;
		else if(allele1 == 0x81) allele1 = 1;
		else allele1 = (allele1 >> 1) + 1;

		if((allele2 >> 1) == 0)  allele2 = 0;
		else if(allele2 == 0x81) allele2 = 1;
		else allele2 = (allele2 >> 1) + 1;

		const U32 internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add, phase);

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

void GenotypeEncoder::updateStatistics(const GenotypeEncoderSlaveHelper& helper){
	if(helper.encoding_type == YON_GT_RLE_DIPLOID_BIALLELIC){
		if(helper.gt_primitive == YON_GT_BYTE){
			++this->stats_.rle_counts[0];
		} else if(helper.gt_primitive == YON_GT_U16){
			++this->stats_.rle_counts[1];
		} else if(helper.gt_primitive == YON_GT_U32){
			++this->stats_.rle_counts[2];
		} else if(helper.gt_primitive == YON_GT_U64){
			++this->stats_.rle_counts[3];
		}
	} else if(helper.encoding_type == YON_GT_RLE_DIPLOID_NALLELIC){
		if(helper.gt_primitive == YON_GT_BYTE){
			++this->stats_.rle_simple_counts[0];
		} else if(helper.gt_primitive == YON_GT_U16){
			++this->stats_.rle_simple_counts[1];
		} else if(helper.gt_primitive == YON_GT_U32){
			++this->stats_.rle_simple_counts[2];
		} else if(helper.gt_primitive == YON_GT_U64){
			++this->stats_.rle_simple_counts[3];
		}
	} else if(helper.encoding_type == YON_GT_BCF_DIPLOID){
		if(helper.gt_primitive == YON_GT_BYTE){
			++this->stats_.diploid_bcf_counts[0];
		} else if(helper.gt_primitive == YON_GT_U16){
			++this->stats_.diploid_bcf_counts[1];
		} else if(helper.gt_primitive == YON_GT_U32){
			++this->stats_.diploid_bcf_counts[2];
		}
	} else if(helper.encoding_type == YON_GT_BCF_STYLE){
		if(helper.gt_primitive == YON_GT_BYTE){
			++this->stats_.bcf_counts[0];
		} else if(helper.gt_primitive == YON_GT_U16){
			++this->stats_.bcf_counts[1];
		} else if(helper.gt_primitive == YON_GT_U32){
			++this->stats_.bcf_counts[2];
		}
	}
}

}
}
