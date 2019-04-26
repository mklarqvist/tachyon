#include "genotype_encoder.h"

namespace tachyon{
namespace algorithm{

GenotypeEncoder::GenotypeEncoder() :
	n_samples(0)
{
}

GenotypeEncoder::GenotypeEncoder(const uint64_t samples) :
	n_samples(samples)
{
}

GenotypeEncoder::~GenotypeEncoder() {}

bool GenotypeEncoder::Encode(const containers::VcfContainer& container,
                             yon1_vnt_t* rcds,
                             block_type& block,
                             const yon_gt_ppa& permutation_array) const
{
	for (uint32_t i = 0; i < container.sizeWithoutCarryOver(); ++i) {
		if (rcds[i].controller.gt_available == false)
			continue;

		GenotypeSummary gt_summary = container.GetGenotypeSummary(i, this->n_samples);

		yon_gt_assess assessed = this->Assess(container[i], gt_summary, permutation_array);
		const uint8_t primitive = assessed.GetCheapestPrimitive();
		//const uint8_t rprimitive = assessed.GetCheapestRawPrimitive();

		// Debug
		/*
		std::cout << rcds[i].rid << "\t" << rcds[i].pos+1 << "\t" << (int)assessed.method << "\t" << (int)primitive << "\t" << assessed.n_cost[primitive] << "\t" << (int)rprimitive-4 << "\t" << assessed.n_cost[rprimitive];
		for (int i = 0; i < 8; ++i) {
			std::cout << "\t" << assessed.n_cost[i];
		}
		for (int i = 0; i < 8; ++i) {
			std::cout << "\t" << assessed.n_runs[i];
		}
		std::cout << std::endl;
		*/

		rcds[i].controller.biallelic         = (container[i]->n_allele == 2);
		rcds[i].controller.diploid           = (gt_summary.base_ploidy == 2);
		rcds[i].controller.gt_has_mixed_phasing = gt_summary.mixed_phasing;
		rcds[i].controller.gt_has_missing    = (gt_summary.n_missing != 0);
		rcds[i].controller.gt_phase_uniform  = gt_summary.phase_if_uniform;
		rcds[i].controller.gt_mixed_ploidy   = (gt_summary.n_vector_end != 0);
		rcds[i].controller.gt_available      = true;
		rcds[i].n_base_ploidy                = gt_summary.base_ploidy;

		if (container[i]->d.fmt[0].n == 2) {
			if (container[i]->n_allele == 2 && gt_summary.n_vector_end == 0) {
				uint64_t n_runs = 0;
				switch(primitive) {
				case(0): n_runs = this->EncodeDiploidBiallelic<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidBiallelic<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT16]); break;
				case(2): n_runs = this->EncodeDiploidBiallelic<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT32]); break;
				case(3): n_runs = this->EncodeDiploidBiallelic<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				rcds[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];

			} else {
				uint64_t n_runs = 0;
				switch(primitive) {
				case(0): n_runs = this->EncodeDiploidMultiAllelic<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidMultiAllelic<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT16]); break;
				case(2): n_runs = this->EncodeDiploidMultiAllelic<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT32]); break;
				case(3): n_runs = this->EncodeDiploidMultiAllelic<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				rcds[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];
			}
		}
		else {
			uint64_t n_runs = 0;
			switch(primitive) {
			case(0): n_runs = this->EncodeMultiploid<uint8_t> (container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT8]);  break;
			case(1): n_runs = this->EncodeMultiploid<uint16_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT16]); break;
			case(2): n_runs = this->EncodeMultiploid<uint32_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT32]); break;
			case(3): n_runs = this->EncodeMultiploid<uint64_t>(container[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT64]); break;
			default:
				std::cerr << "illegal primitive type" << std::endl;
				return false;
			}

			rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
			rcds[i].controller.gt_compression_type = YON_GT_RLE_NPLOID;
			block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
			++block.base_containers[YON_BLK_GT_SUPPORT];
		}
	}
	return true;
}

bool GenotypeEncoder::Encode(yon1_vnt_t* rcds,
                             const uint32_t n_rcds,
                             block_type& block,
                             const yon_gt_ppa& permutation_array) const
{
	for (uint32_t i = 0; i < n_rcds; ++i) {
		if (rcds[i].controller.gt_available == false)
			continue;

		GenotypeSummary gt_summary;
		gt_summary.Evaluate(rcds[i], this->n_samples);

		yon_gt_assess assessed = this->Assess(rcds[i], gt_summary, permutation_array);
		const uint8_t primitive = assessed.GetCheapestPrimitive();
		//const uint8_t rprimitive = assessed.GetCheapestRawPrimitive();

		// Debug
		/*
		std::cout << rcds[i].rid << "\t" << rcds[i].pos+1 << "\t" << (int)assessed.method << "\t" << (int)primitive << "\t" << assessed.n_cost[primitive] << "\t" << (int)rprimitive-4 << "\t" << assessed.n_cost[rprimitive];
		for (int i = 0; i < 8; ++i) {
			std::cout << "\t" << assessed.n_cost[i];
		}
		for (int i = 0; i < 8; ++i) {
			std::cout << "\t" << assessed.n_runs[i];
		}
		std::cout << std::endl;
		*/

		rcds[i].controller.biallelic         = (rcds[i].n_alleles == 2);
		rcds[i].controller.diploid           = (gt_summary.base_ploidy == 2);
		rcds[i].controller.gt_has_mixed_phasing = gt_summary.mixed_phasing;
		rcds[i].controller.gt_has_missing    = (gt_summary.n_missing != 0);
		rcds[i].controller.gt_phase_uniform  = gt_summary.phase_if_uniform;
		rcds[i].controller.gt_mixed_ploidy   = (gt_summary.n_vector_end != 0);
		rcds[i].controller.gt_available      = true;
		rcds[i].n_base_ploidy                = gt_summary.base_ploidy;

		if (rcds[i].gt->m == 2) {
			if (rcds[i].n_alleles == 2 && gt_summary.n_vector_end == 0) {
				uint64_t n_runs = 0;
				switch(primitive) {
				case(0): n_runs = this->EncodeDiploidBiallelic<uint8_t> (rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidBiallelic<uint16_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT16]); break;
				case(2): n_runs = this->EncodeDiploidBiallelic<uint32_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT32]); break;
				case(3): n_runs = this->EncodeDiploidBiallelic<uint64_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				rcds[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_BIALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];

			} else {
				uint64_t n_runs = 0;
				switch(primitive) {
				case(0): n_runs = this->EncodeDiploidMultiAllelic<uint8_t> (rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT8]);  break;
				case(1): n_runs = this->EncodeDiploidMultiAllelic<uint16_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT16]); break;
				case(2): n_runs = this->EncodeDiploidMultiAllelic<uint32_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT32]); break;
				case(3): n_runs = this->EncodeDiploidMultiAllelic<uint64_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_S_INT64]); break;
				default:
					std::cerr << "illegal primitive type" << std::endl;
					return false;
				}

				rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
				rcds[i].controller.gt_compression_type = YON_GT_RLE_DIPLOID_NALLELIC;
				block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
				++block.base_containers[YON_BLK_GT_SUPPORT];
			}
		}
		else {
			uint64_t n_runs = 0;
			switch(primitive) {
			case(0): n_runs = this->EncodeMultiploid<uint8_t> (rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT8]);  break;
			case(1): n_runs = this->EncodeMultiploid<uint16_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT16]); break;
			case(2): n_runs = this->EncodeMultiploid<uint32_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT32]); break;
			case(3): n_runs = this->EncodeMultiploid<uint64_t>(rcds[i],gt_summary,permutation_array,block.base_containers[YON_BLK_GT_N_INT64]); break;
			default:
				std::cerr << "illegal primitive type" << std::endl;
				return false;
			}

			rcds[i].controller.gt_primtive_type    = TACHYON_GT_PRIMITIVE_TYPE(primitive);
			rcds[i].controller.gt_compression_type = YON_GT_RLE_NPLOID;
			block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)n_runs);
			++block.base_containers[YON_BLK_GT_SUPPORT];
		}
	}
	return true;
}

yon_gt_assess GenotypeEncoder::Assess(const bcf1_t* entry,
                                      const GenotypeSummary& gt_summary,
                                      const yon_gt_ppa& permutation_array) const
{
	// Special case of diploid record.
	yon_gt_assess assessed;
	if (entry->d.fmt[0].n == 2) {
		if (entry->n_allele == 2 && gt_summary.n_vector_end == 0)
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

yon_gt_assess GenotypeEncoder::Assess(const yon1_vnt_t& entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& permutation_array) const {
	if (entry.controller.gt_available == false)
		return yon_gt_assess();

	if (entry.gt == nullptr) {
		std::cerr << "no gt set" << std::endl;
	}

	if (entry.gt->d_exp == nullptr) {
		std::cerr << "data not expanded" << std::endl;
	}

	yon_gt_assess assessed;
	if (entry.gt->m == 2) {
		if (entry.n_alleles == 2)
			assessed = this->AssessDiploidBiallelic(entry, gt_summary, permutation_array);
		else
			assessed = this->AssessDiploidMultiAllelic(entry, gt_summary, permutation_array);
	} else {
		assessed = this->AssessMultiploid(entry, gt_summary, permutation_array);
	}

	return assessed;
}

yon_gt_assess GenotypeEncoder::AssessDiploidBiallelic(const bcf1_t* entry,
                                                      const GenotypeSummary& gt_summary,
                                                      const yon_gt_ppa& permutation_array) const
{
	assert(entry->d.fmt[0].n == 2);
	assert(entry->n_allele == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	// Track all possible outcomes.
	// 1: uint8_t + Permuted
	// 2: uint16_t  + Permuted
	// 3: uint32_t  + Permuted
	// 4: uint64_t  + Permuted
	// 5: uint8_t + No permutation
	// 6: uint16_t  + No permutation
	// 7: uint32_t  + No permutation
	// 8: uint64_t  + No permutation
	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;

	// 1 + hasMissing + hasMixedPhasing
	const uint8_t shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const uint8_t add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	uint64_t limits[4];
	limits[0] = pow(2, 8*sizeof(uint8_t) - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(uint16_t)  - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(uint32_t)  - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(uint64_t)  - (base_ploidy*shift + add)) - 1;

	uint32_t rle_current_ref     = YON_PACK_GT_DIPLOID(gt[0], gt[1], shift, add);
	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID(gt[permutation_array[0] * sizeof(int8_t) * base_ploidy],
	                                              gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1],
												  shift, add);

	// Keep track of the linear offset in the genotype
	// data stream. The permuted offset is computed directly
	// and does not need to be tracked.
	uint32_t l_gt_offset = 2;

	// Iterate over all available samples.
	for (uint32_t i = 1; i < this->n_samples; ++i, l_gt_offset += 2) {
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		uint32_t rle_current     = YON_PACK_GT_DIPLOID(gt[l_gt_offset], gt[l_gt_offset + 1], shift, add);
		uint32_t rle_ppa_current = YON_PACK_GT_DIPLOID(gt_ppa_target[0], gt_ppa_target[1], shift, add);

		if (rle_current != rle_current_ref) {
			for (uint32_t k = 4; k < 8; ++k) ++n_runs[k];
			for (uint32_t k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 4; k < 8; ++k) {
			if (l_runs[k] == limits[k-4]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if (rle_ppa_current != rle_ppa_current_ref) {
			for (uint32_t k = 0; k < 4; ++k) ++n_runs[k];
			for (uint32_t k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a breakAssessDiploidBiallelic
		for (uint32_t k = 0; k < 4; ++k) {
			if (l_runs[k] == limits[k]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

	}
	assert(l_gt_offset == l_gt);
	for (uint32_t k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for (uint32_t k = 0; k < 4; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*(k+1);
	}
	for (uint32_t k = 4; k < 8; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*((k-4)+1);
	}
	sum.method = 0;

	/*
	std::cout << entry->pos + 1 << "\tR";
	for (uint32_t i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for (uint32_t i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessDiploidBiallelic(const yon1_vnt_t& entry,
                                                      const GenotypeSummary& gt_summary,
                                                      const yon_gt_ppa& ppa) const
{
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	assert(entry.gt->m == 2);
	assert(entry.n_alleles == 2);
	assert(ppa.n_s == entry.gt->n_s);
	const uint8_t base_ploidy = entry.gt->m;

	// Track all possible outcomes.
	// 1: uint8_t + Permuted
	// 2: uint16_t  + Permuted
	// 3: uint32_t  + Permuted
	// 4: uint64_t  + Permuted
	// 5: uint8_t + No permutation
	// 6: uint16_t  + No permutation
	// 7: uint32_t  + No permutation
	// 8: uint64_t  + No permutation
	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;

	// 1 + hasMissing + hasMixedPhasing
	const uint8_t shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const uint8_t add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	uint64_t limits[4];
	limits[0] = pow(2, 8*sizeof(uint8_t)  - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(uint16_t) - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(uint32_t) - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(uint64_t) - (base_ploidy*shift + add)) - 1;

	uint32_t rle_current_ref     = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[0], shift, add);
	uint32_t rle_ppa_current_ref = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[ppa[0]], shift, add);

	// Iterate over all available samples.
	for (uint32_t i = 1; i < this->n_samples; ++i) {
		uint32_t rle_current     = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[i], shift, add);
		uint32_t rle_ppa_current = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[ppa[i]], shift, add);

		if (rle_current != rle_current_ref) {
			for (uint32_t k = 4; k < 8; ++k) ++n_runs[k];
			for (uint32_t k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 4; k < 8; ++k) {
			if (l_runs[k] == limits[k-4]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if (rle_ppa_current != rle_ppa_current_ref) {
			for (uint32_t k = 0; k < 4; ++k) ++n_runs[k];
			for (uint32_t k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a breakAssessDiploidBiallelic
		for (uint32_t k = 0; k < 4; ++k) {
			if (l_runs[k] == limits[k]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

	}

	for (uint32_t k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for (uint32_t k = 0; k < 4; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*(k+1);
	}
	for (uint32_t k = 4; k < 8; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*((k-4)+1);
	}
	sum.method = 0;

	/*
	std::cout << entry.rid << ":" << entry.pos + 1;
	for (uint32_t i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1);
	for (uint32_t i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1);
	std::cout << std::endl;
	*/

	return yon_gt_assess();
}

yon_gt_assess GenotypeEncoder::AssessDiploidMultiAllelic(const bcf1_t* entry,
                                                         const GenotypeSummary& gt_summary,
                                                         const yon_gt_ppa& permutation_array) const
{
	assert(entry->d.fmt[0].n == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	// Track all possible outcomes.
	// 1: uint8_t + Permuted
	// 2: uint16_t  + Permuted
	// 3: uint32_t  + Permuted
	// 4: uint64_t  + Permuted
	// 5: uint8_t + No permutation
	// 6: uint16_t  + No permutation
	// 7: uint32_t  + No permutation
	// 8: uint64_t  + No permutation
	int64_t n_runs[8]; // Number of runs.
	int64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;

	// Assess RLE cost
	const uint8_t shift = ceil(log2(entry->n_allele + 2 + 1));
	const uint8_t add   = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size.
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
	int64_t limits[4];
	limits[0] = pow(2, 8*sizeof(uint8_t) - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(uint16_t)  - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(uint32_t)  - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(uint64_t)  - (base_ploidy*shift + add)) - 1;
	bool banned_limit[4];
	memset(banned_limit, 0, 4*sizeof(bool));
	if (limits[0] <= 0) { limits[0] = std::numeric_limits<int64_t>::max(); banned_limit[0] = true; }
	if (limits[1] <= 0) { limits[1] = std::numeric_limits<int64_t>::max(); banned_limit[1] = true; }
	if (limits[2] <= 0) { limits[2] = std::numeric_limits<int64_t>::max(); banned_limit[2] = true; }

	uint8_t gt_remap[256];
	memset(gt_remap, 256, 255);
	for (uint32_t i = 0; i <= entry->n_allele; ++i) {
		gt_remap[i << 1]       = ((i+1) << 1);
		gt_remap[(i << 1) + 1] = ((i+1) << 1) + 1;
	}
	gt_remap[0]   = 0;
	gt_remap[129] = 1;

	uint32_t rle_current_ref     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[0]] >> 1,
	                                                       gt_remap[gt[1]] >> 1,
	                                                       shift, add,
	                                                       gt_remap[gt[1]]);

	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy]] >> 1,
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
	for (uint32_t i = 1; i < this->n_samples; ++i, l_gt_offset += base_ploidy) {
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		uint32_t rle_current     = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[l_gt_offset]] >> 1, gt_remap[gt[l_gt_offset+1]] >> 1, shift, add, gt_remap[gt[1]]);
		uint32_t rle_ppa_current = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt_ppa_target[0]] >> 1,
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

		if (rle_current != rle_current_ref) {
			for (uint32_t k = 4; k < 8; ++k) ++n_runs[k];
			for (uint32_t k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 4; k < 8; ++k) {
			if (l_runs[k] == limits[k-4]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if (rle_ppa_current != rle_ppa_current_ref) {
			for (uint32_t k = 0; k < 4; ++k) ++n_runs[k];
			for (uint32_t k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 0; k < 4; ++k) {
			if (l_runs[k] == limits[k]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}
	}
	assert(l_gt_offset == l_gt);
	for (uint32_t k = 0; k < 8; ++k) ++n_runs[k];

	uint8_t banned_sum = 0;
	for (int i = 0; i < 4; ++i) banned_sum += banned_limit[i];
	assert(banned_sum != 4);

	yon_gt_assess sum;
	for (uint32_t k = 0; k < 4; ++k) {
		if (banned_limit[k]) {
			sum.n_runs[k] = std::numeric_limits<uint64_t>::max();
			sum.n_cost[k] = std::numeric_limits<uint64_t>::max();
		} else {
			sum.n_runs[k] = n_runs[k];
			sum.n_cost[k] = n_runs[k]*(k+1);
		}
	}
	for (uint32_t k = 4; k < 8; ++k) {
		if (banned_limit[k-4]) {
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
	for (uint32_t i = 0; i < 8; ++i)
		std::cout << "\t" << sum.n_runs[i] << "\t" << sum.n_cost[i];
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessDiploidMultiAllelic(const yon1_vnt_t& entry,
                                                         const GenotypeSummary& gt_summary,
                                                         const yon_gt_ppa& ppa) const
{
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	assert(entry.gt->m == 2);
	assert(ppa.n_s == entry.gt->n_s);
	assert(this->n_samples == entry.gt->n_s);
	const uint8_t base_ploidy = entry.gt->m;

	// Track all possible outcomes.
	// 1: uint8_t + Permuted
	// 2: uint16_t  + Permuted
	// 3: uint32_t  + Permuted
	// 4: uint64_t  + Permuted
	// 5: uint8_t + No permutation
	// 6: uint16_t  + No permutation
	// 7: uint32_t  + No permutation
	// 8: uint64_t  + No permutation
	int64_t n_runs[8]; // Number of runs.
	int64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;

	// Assess RLE cost
	const uint8_t shift = ceil(log2(entry.n_alleles + 2 + 1));
	const uint8_t add   = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size.
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
	int64_t limits[4];
	limits[0] = pow(2, 8*sizeof(uint8_t)  - (base_ploidy*shift + add)) - 1;
	limits[1] = pow(2, 8*sizeof(uint16_t) - (base_ploidy*shift + add)) - 1;
	limits[2] = pow(2, 8*sizeof(uint32_t) - (base_ploidy*shift + add)) - 1;
	limits[3] = pow(2, 8*sizeof(uint64_t) - (base_ploidy*shift + add)) - 1;
	bool banned_limit[4];
	memset(banned_limit, 0, 4*sizeof(bool));
	if (limits[0] <= 0) { limits[0] = std::numeric_limits<int64_t>::max(); banned_limit[0] = true; }
	if (limits[1] <= 0) { limits[1] = std::numeric_limits<int64_t>::max(); banned_limit[1] = true; }
	if (limits[2] <= 0) { limits[2] = std::numeric_limits<int64_t>::max(); banned_limit[2] = true; }

	uint32_t rle_current_ref = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[0], shift, add);
	uint32_t rle_ppa_current_ref = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[ppa[0]], shift, add);

	assert(ceil(log2(entry.gt->d_exp[0].allele[0] >> 1)) <= shift);
	assert(ceil(log2(entry.gt->d_exp[0].allele[1] >> 1)) <= shift);
	assert(ceil(log2(entry.gt->d_exp[ppa[0]].allele[0] >> 1)) <= shift);
	assert(ceil(log2(entry.gt->d_exp[ppa[0]].allele[0] >> 1)) <= shift);

	// Iterate over all available samples.
	for (uint32_t i = 1; i < this->n_samples; ++i) {
		uint32_t rle_current     = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[i], shift, add);
		uint32_t rle_ppa_current = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[ppa[i]], shift, add);

		assert(ceil(log2(entry.gt->d_exp[i].allele[0] >> 1)) <= shift);
		assert(ceil(log2(entry.gt->d_exp[i].allele[1] >> 1)) <= shift);
		assert(ceil(log2(entry.gt->d_exp[ppa[i]].allele[0] >> 1)) <= shift);
		assert(ceil(log2(entry.gt->d_exp[ppa[i]].allele[1] >> 1)) <= shift);

		if (rle_current != rle_current_ref) {
			for (uint32_t k = 4; k < 8; ++k) ++n_runs[k];
			for (uint32_t k = 4; k < 8; ++k) l_runs[k] = 0;
			rle_current_ref = rle_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 4; k < 8; ++k) {
			if (l_runs[k] == limits[k-4]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if (rle_ppa_current != rle_ppa_current_ref) {
			for (uint32_t k = 0; k < 4; ++k) ++n_runs[k];
			for (uint32_t k = 0; k < 4; ++k) l_runs[k] = 0;
			rle_ppa_current_ref = rle_ppa_current;
		}

		// Overflow: trigger a break
		for (uint32_t k = 0; k < 4; ++k) {
			if (l_runs[k] == limits[k]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}
	}

	for (uint32_t k = 0; k < 8; ++k) ++n_runs[k];

	uint8_t banned_sum = 0;
	for (int i = 0; i < 4; ++i) banned_sum += banned_limit[i];
	assert(banned_sum != 4);

	yon_gt_assess sum;
	for (uint32_t k = 0; k < 4; ++k) {
		if (banned_limit[k]) {
			sum.n_runs[k] = std::numeric_limits<uint64_t>::max();
			sum.n_cost[k] = std::numeric_limits<uint64_t>::max();
		} else {
			sum.n_runs[k] = n_runs[k];
			sum.n_cost[k] = n_runs[k]*(k+1);
		}
	}
	for (uint32_t k = 4; k < 8; ++k) {
		if (banned_limit[k-4]) {
			sum.n_runs[k] = std::numeric_limits<uint64_t>::max();
			sum.n_cost[k] = std::numeric_limits<uint64_t>::max();
		} else {
			sum.n_runs[k] = n_runs[k];
			sum.n_cost[k] = n_runs[k]*((k-4)+1);
		}
	}
	sum.method = 1;

	/*
	std::cout << entry.rid << ":" << entry.pos + 1;
	for (uint32_t i = 0; i < 8; ++i)
		std::cout << "\t" << sum.n_runs[i] << "\t" << sum.n_cost[i];
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessMultiploid(const bcf1_t* entry,
                                                const GenotypeSummary& gt_summary,
                                                const yon_gt_ppa& permutation_array) const
{
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;
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
	for (uint32_t i = 1; i < this->n_samples; ++i, l_gt_offset += base_ploidy) {
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
		uint64_t hash_value     = XXH64(&gt[l_gt_offset], sizeof(int8_t) * base_ploidy, 89231478);
		uint64_t hash_value_ppa = XXH64(gt_ppa_target,    sizeof(int8_t) * base_ploidy, 89231478);

		if (hash_value_ppa != hash_value_ppa_ref) {
			for (uint32_t k = 0; k < 4; ++k) ++n_runs[k];
			for (uint32_t k = 0; k < 4; ++k) l_runs[k] = 0;
			hash_value_ppa_ref = hash_value_ppa;
		}

		// Overflow: trigger a break
		for (uint32_t k = 0; k < 4; ++k) {
			if (l_runs[k] == limits[k]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}

		if (hash_value != hash_value_ref) {
			for (uint32_t k = 4; k < 8; ++k) ++n_runs[k];
			for (uint32_t k = 4; k < 8; ++k) l_runs[k] = 0;
			hash_value_ref = hash_value;
		}

		// Overflow: trigger a break
		for (uint32_t k = 4; k < 8; ++k) {
			if (l_runs[k] == limits[k-4]) { ++n_runs[k]; l_runs[k] = 0; }
			++l_runs[k];
		}
	}
	assert(l_gt_offset == l_gt);
	for (uint32_t k = 0; k < 8; ++k) ++n_runs[k];

	yon_gt_assess sum;
	for (uint32_t k = 0; k < 4; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*(k+1);
	}
	for (uint32_t k = 4; k < 8; ++k) {
		sum.n_runs[k] = n_runs[k];
		sum.n_cost[k] = n_runs[k]*((k-4)+1);
	}
	sum.method = 2;

	/*
	std::cout << entry->pos + 1 << "\tX\t" << this->n_samples;
	for (uint32_t i = 0; i < 4; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*(i+1)*base_ploidy;
	for (uint32_t i = 4; i < 8; ++i)
		std::cout << "\t" << n_runs[i] << "\t" << n_runs[i]*((i-4)+1)*base_ploidy;
	std::cout << std::endl;
	*/

	return sum;
}

yon_gt_assess GenotypeEncoder::AssessMultiploid(const yon1_vnt_t& entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const {
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	const uint8_t   base_ploidy = entry.gt->m;
	assert(ppa.n_s == entry.gt->n_s);

	uint64_t n_runs[8]; // Number of runs.
	uint64_t l_runs[8]; // Current run length.
	for (uint32_t i = 0; i < 8; ++i) l_runs[i] = 1;
	for (uint32_t i = 0; i < 8; ++i) n_runs[i] = 0;
	uint64_t limits[4];
	limits[0] = std::numeric_limits<uint8_t>::max();
	limits[1] = std::numeric_limits<uint16_t>::max();
	limits[2] = std::numeric_limits<uint32_t>::max();
	limits[3] = std::numeric_limits<uint64_t>::max();

	std::cerr << "not implemented yet" << std::endl;
	exit(1);

	return yon_gt_assess();
}

}
}
