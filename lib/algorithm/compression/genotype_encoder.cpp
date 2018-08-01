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

bool GenotypeEncoder::Encode(const containers::VcfContainer& container, meta_type* meta_entries, block_type& block, const algorithm::yon_gt_ppa& permutation_array) const{
	for(U32 i = 0; i < container.sizeWithoutCarryOver(); ++i){
		io::VcfGenotypeSummary gt_summary = container.GetGenotypeSummary(i, this->n_samples);

		yon_gt_assess assessed = this->Assess(container[i], gt_summary, permutation_array);
		const uint8_t primitive = assessed.GetCheapestPrimitive();

		meta_entries[i].controller.biallelic        = (container[i]->n_allele == 2);
		meta_entries[i].controller.diploid          = (gt_summary.base_ploidy == 2);
		meta_entries[i].controller.gt_mixed_phasing = gt_summary.mixed_phasing;
		meta_entries[i].controller.gt_anyMissing    = (gt_summary.n_missing != 0);
		meta_entries[i].controller.gt_anyNA         = (gt_summary.n_vector_end != 0);
		meta_entries[i].controller.gt_phase         = gt_summary.phase_if_uniform;
		meta_entries[i].controller.mixed_ploidy     = (gt_summary.n_vector_end != 0);
		meta_entries[i].controller.gt_available     = true;

		if(container[i]->d.fmt[0].n == 2){
			if(container[i]->n_allele == 2 && gt_summary.n_vector_end == 0){
				uint64_t n_runs = 0;
				switch(primitive){
				case(0): n_runs = this->EncodeDiploidBiallelic<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidBiallelic<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT16]); break;
				case(2): n_runs = this->EncodeDiploidBiallelic<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT32]); break;
				case(3): n_runs = this->EncodeDiploidBiallelic<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				meta_entries[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				meta_entries[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];

			} else {
				uint64_t n_runs = 0;
				switch(primitive){
				case(0): n_runs = this->EncodeDiploidMultiAllelic<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidMultiAllelic<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT16]); break;
				case(2): n_runs = this->EncodeDiploidMultiAllelic<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT32]); break;
				case(3): n_runs = this->EncodeDiploidMultiAllelic<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				meta_entries[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				meta_entries[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];
			}
		}
		else {
			uint64_t n_runs = 0;
			switch(primitive){
			case(0): n_runs = this->EncodeMultiploid<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT8]);  break;
			case(1): n_runs = this->EncodeMultiploid<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT16]); break;
			case(2): n_runs = this->EncodeMultiploid<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT32]); break;
			case(3): n_runs = this->EncodeMultiploid<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT64]); break;
			default:
				std::cerr << "illegal primitive type" << std::endl;
				return false;
			}

			meta_entries[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
			meta_entries[i].controller.gt_compression_type = YON_GT_BCF_STYLE;
			block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)n_runs);
			++block.base_containers[YON_BLK_GT_SUPPORT];
		}
	}
	return true;
}

yon_gt_assess GenotypeEncoder::Assess(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
	// Special case of diploid record.
	yon_gt_assess assessed;
	if(entry->d.fmt[0].n == 2){
		if(entry->n_allele == 2 && gt_summary.n_vector_end == 0)
			assessed = this->AssessDiploidBiallelic(entry,gt_summary,permutation_array);
		else
			assessed = this->AssessDiploidMultiAllelic(entry,gt_summary,permutation_array);
	}
	// All other ploidy is assessed with this multiploid function.
	else {
		assessed = this->AssessMultiploid(entry,gt_summary,permutation_array);
	}
	return assessed;
}

yon_gt_assess GenotypeEncoder::AssessDiploidBiallelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
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
	for(U32 i = 0; i < 8; ++i) l_runs[i] = 1;
	for(U32 i = 0; i < 8; ++i) n_runs[i] = 0;

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
			for(U32 k = 4; k < 8; ++k) ++n_runs[k];
			for(U32 k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for(U32 k = 4; k < 8; ++k){
			if(l_runs[k] == limits[k-4]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if(rle_ppa_current != rle_ppa_current_ref){
			for(U32 k = 0; k < 4; ++k) ++n_runs[k];
			for(U32 k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for(U32 k = 0; k < 4; ++k){
			if(l_runs[k] == limits[k]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

	}
	assert(l_gt_offset == l_gt);
	for(U32 k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for(U32 k = 0; k < 4; ++k){
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*(k+1);
	}
	for(U32 k = 4; k < 8; ++k){
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*((k-4)+1);
	}
	sum.method = 0;

	/*
	std::cout << entry->pos + 1 << "\tR";
	for(U32 i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for(U32 i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessDiploidMultiAllelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
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
	for(U32 i = 0; i < 8; ++i) l_runs[i] = 1;
	for(U32 i = 0; i < 8; ++i) n_runs[i] = 0;

	// Assess RLE cost
	const BYTE shift = ceil(log2(entry->n_allele + 2 + 1));
	const BYTE add   = gt_summary.mixed_phasing  ? 1 : 0;

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
	bool banned_limit[4];
	if(limits[0] <= 0){ limits[0] = std::numeric_limits<int64_t>::max(); banned_limit[0] = true; }
	if(limits[1] <= 0){ limits[1] = std::numeric_limits<int64_t>::max(); banned_limit[1] = true; }
	if(limits[2] <= 0){ limits[2] = std::numeric_limits<int64_t>::max(); banned_limit[2] = true; }

	uint8_t gt_remap[256];
	memset(gt_remap, 256, 255);
	for(U32 i = 0; i <= entry->n_allele; ++i){
		gt_remap[i << 1]       = (i << 1) + 1;
		gt_remap[(i << 1) + 1] = (i << 1) + 2;
	}
	gt_remap[0] = 0;
	gt_remap[129] = 1;

	U32 rle_current_ref     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[0]] >> 1, gt_remap[gt[1]] >> 1, shift, add, gt_remap[gt[1]]);
	U32 rle_ppa_current_ref = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy]] >> 1,
												           gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]] >> 1,
												           shift, add,
												           gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]]);

	assert(ceil(log2(gt_remap[gt[0]] >> 1)) <= shift);
	assert(ceil(log2(gt_remap[gt[1]] >> 1)) <= shift);
	assert(ceil(log2(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy]] >> 1)) <= shift);
	assert(ceil(log2(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]] >> 1)) <= shift);

	// Keep track of the linear offset in the genotype
	// data stream. The permuted offset is computed directly
	// and does not need to be tracked.
	uint32_t l_gt_offset = base_ploidy;

	//std::cerr << entry->pos + 1 << "\t" << (gt[0]>>1) << "|" << (gt[1]>>1);

	// Iterate over all available samples.
	for(U32 i = 1; i < this->n_samples; ++i, l_gt_offset += base_ploidy){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		U32 rle_current     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[l_gt_offset]] >> 1, gt_remap[gt[l_gt_offset+1]] >> 1, shift, add, gt_remap[gt[1]]);
		U32 rle_ppa_current = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt_ppa_target[0]] >> 1,
		                                                   gt_remap[gt_ppa_target[1]] >> 1,
												           shift, add,
														   gt_remap[gt_ppa_target[1]]);

		assert(ceil(log2(gt_remap[gt[l_gt_offset]] >> 1)) <= shift);
		assert(ceil(log2(gt_remap[gt[l_gt_offset+1]] >> 1)) <= shift);
		assert(ceil(log2(gt_remap[gt_ppa_target[0]] >> 1)) <= shift);
		assert(ceil(log2(gt_remap[gt_ppa_target[1]] >> 1)) <= shift);

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
			for(U32 k = 4; k < 8; ++k) ++n_runs[k];
			for(U32 k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for(U32 k = 4; k < 8; ++k){
			if(l_runs[k] == limits[k-4]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if(rle_ppa_current != rle_ppa_current_ref){
			for(U32 k = 0; k < 4; ++k) ++n_runs[k];
			for(U32 k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for(U32 k = 0; k < 4; ++k){
			if(l_runs[k] == limits[k]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}
	}
	//std::cerr << std::endl;
	assert(l_gt_offset == l_gt);
	for(U32 k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for(U32 k = 0; k < 4; ++k){
		if(banned_limit[k]){
			sum.n_runs[k] = std::numeric_limits<uint64_t>::max();
			sum.n_cost[k] = std::numeric_limits<uint64_t>::max();
		} else {
			sum.n_runs[k] = n_runs[k];
			sum.n_cost[k] = n_runs[k]*(k+1);
		}
	}
	for(U32 k = 4; k < 8; ++k){
		if(banned_limit[k]){
			sum.n_runs[k] = std::numeric_limits<uint64_t>::max();
			sum.n_cost[k] = std::numeric_limits<uint64_t>::max();
		} else {
			sum.n_runs[k] = n_runs[k];
			sum.n_cost[k] = n_runs[k]*((k-4)+1);
		}
	}
	sum.method = 1;

	/*
	std::cout << entry->pos + 1 << "\tM";
	for(U32 i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for(U32 i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessMultiploid(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const{
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_samples * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	for(U32 i = 0; i < 8; ++i) l_runs[i] = 1;
	for(U32 i = 0; i < 8; ++i) n_runs[i] = 0;
	uint64_t limits[4];
	limits[0] = std::numeric_limits<uint8_t>::max();
	limits[1] = std::numeric_limits<uint16_t>::max();
	limits[2] = std::numeric_limits<uint32_t>::max();
	limits[3] = std::numeric_limits<uint64_t>::max();

	//std::cerr << entry->pos + 1 << " : " << this->n_samples << "\t";
	uint64_t hash_value_ref = XXH64(&gt[0], sizeof(int8_t) * base_ploidy, 89231478);
	uint64_t hash_value_ppa_ref = XXH64(&gt[permutation_array[0] * sizeof(int8_t) * base_ploidy], sizeof(int8_t) * base_ploidy, 89231478);

	uint32_t l_gt_offset = base_ploidy;

	// Iterate over all available samples.
	for(U32 i = 1; i < this->n_samples; ++i, l_gt_offset += base_ploidy){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		uint64_t hash_value     = XXH64(&gt[l_gt_offset], sizeof(int8_t) * base_ploidy, 89231478);
		uint64_t hash_value_ppa = XXH64(gt_ppa_target,    sizeof(int8_t) * base_ploidy, 89231478);

		if(hash_value_ppa != hash_value_ppa_ref){
			for(U32 k = 0; k < 4; ++k) ++n_runs[k];
			for(U32 k = 0; k < 4; ++k) l_runs[k] = 0;
			hash_value_ppa_ref = hash_value_ppa;
		}

		// Overflow: trigger a break
		for(U32 k = 0; k < 4; ++k){
			if(l_runs[k] == limits[k]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if(hash_value != hash_value_ref){
			for(U32 k = 4; k < 8; ++k) ++n_runs[k];
			for(U32 k = 4; k < 8; ++k) l_runs[k] = 0;
			hash_value_ref = hash_value;
		}

		// Overflow: trigger a break
		for(U32 k = 4; k < 8; ++k){
			if(l_runs[k] == limits[k-4]){ ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}
	}
	assert(l_gt_offset == l_gt);
	for(U32 k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for(U32 k = 0; k < 4; ++k){
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*(k+1);
	}
	for(U32 k = 4; k < 8; ++k){
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*((k-4)+1);
	}
	sum.method = 2;

	/*
	std::cout << entry->pos + 1 << "\tX\t" << this->n_samples;
	for(U32 i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1)*base_ploidy;
	for(U32 i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1)*base_ploidy;
	std::cout << std::endl;
	*/

	return sum;
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
