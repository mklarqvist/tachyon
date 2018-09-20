#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>
#include <thread>

#include "containers/variant_containers.h"
#include "core/genotypes.h"
#include "io/vcf_utils.h"
#include "core/variant_record.h"
#include "containers/vcf_container.h"

#include "third_party/xxhash/xxhash.h"

namespace tachyon{
namespace algorithm{

const uint8_t BCF_UNPACK_TACHYON[3] = {2, 0, 1};
const uint8_t YON_RCD_DIPLOID_BIALLELIC[4] = {2, 0, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TACHYON[((A) >> 1)]
#define YON_PACK_GT_DIPLOID(A, B, SHIFT, ADD) (BCF_UNPACK_GENOTYPE(A) << ((SHIFT) + (ADD))) | (BCF_UNPACK_GENOTYPE(B) << (ADD)) | ((A) & (ADD))
#define YON_PACK_GT_DIPLOID_NALLELIC(A, B, SHIFT, ADD, PHASE) ((A) << ((SHIFT) + (ADD))) | ((B) << (ADD)) | ((PHASE) & (ADD))
#define YON_PACK_GT_RCD_DIPLOID(A, B, SHIFT, ADD) (YON_RCD_DIPLOID_BIALLELIC[YON_GT_RCD_ALLELE_UNPACK(A)] << ((SHIFT) + (ADD))) | (YON_RCD_DIPLOID_BIALLELIC[YON_GT_RCD_ALLELE_UNPACK(B)] << (ADD)) | ((A) & (ADD))
#define YON_PACK_GT_RCD_DIPLOID_EXPAND(A, SHIFT, ADD) (YON_RCD_DIPLOID_BIALLELIC[YON_GT_RCD_ALLELE_UNPACK(A.allele[0])] << ((SHIFT) + (ADD))) | (YON_RCD_DIPLOID_BIALLELIC[YON_GT_RCD_ALLELE_UNPACK(A.allele[1])] << (ADD)) | ((A.allele[1]) & (ADD))
#define YON_PACK_GT_RCD_NALLELIC(A, B, SHIFT, ADD, PHASE) (YON_GT_RCD_ALLELE_UNPACK(A) << ((SHIFT) + (ADD))) | (YON_GT_RCD_ALLELE_UNPACK(B) << (ADD)) | ((PHASE) & (ADD))
#define YON_PACK_GT_RCD_NALLELIC_EXPAND(A, SHIFT, ADD) (YON_GT_RCD_ALLELE_UNPACK(A.allele[0]) << ((SHIFT) + (ADD))) | (YON_GT_RCD_ALLELE_UNPACK(A.allele[1]) << (ADD)) | ((A.allele[1]) & (ADD))

struct yon_gt_assess {
	uint8_t GetCheapestPrimitive(void) const {
		uint64_t n_cost_best = this->n_cost[0];
		uint8_t best_primitive = 0;
		for(int i = 1; i < 4; ++i){
			if(this->n_cost[i] < n_cost_best){
				n_cost_best = this->n_cost[i];
				best_primitive = i;
			}
		}
		return(best_primitive);
	}

	uint8_t GetCheapestRawPrimitive(void) const {
		uint64_t n_cost_best = this->n_cost[4];
		uint8_t best_primitive = 4;
		for(int i = 5; i < 8; ++i){
			if(this->n_cost[i] < n_cost_best){
				n_cost_best = this->n_cost[i];
				best_primitive = i;
			}
		}
		return(best_primitive);
	}


	uint8_t  method;
	uint64_t n_runs[8];
	uint64_t n_cost[8];
};

struct GenotypeEncoderStatistics{
	GenotypeEncoderStatistics(){
		memset(this->rle_counts,         0, sizeof(uint64_t)*4);
		memset(this->rle_simple_counts,  0, sizeof(uint64_t)*4);
		memset(this->diploid_bcf_counts, 0, sizeof(uint64_t)*3);
		memset(this->bcf_counts,         0, sizeof(uint64_t)*3);
	}

	uint64_t getTotal(void) const{
		uint64_t total = 0;
		for(uint32_t i = 0; i < 4; ++i) total += this->rle_counts[i];
		for(uint32_t i = 0; i < 4; ++i) total += this->rle_simple_counts[i];
		for(uint32_t i = 0; i < 3; ++i) total += this->diploid_bcf_counts[i];
		for(uint32_t i = 0; i < 3; ++i) total += this->bcf_counts[i];
		return(total);
	}

	uint64_t rle_counts[4];
	uint64_t rle_simple_counts[4];
	uint64_t diploid_bcf_counts[3];
	uint64_t bcf_counts[3];
};

class GenotypeEncoder {
public:
	typedef GenotypeEncoder              self_type;
	typedef io::BasicBuffer              buffer_type;
	typedef yon1_dc_t    container_type;
	typedef yon1_vb_t     block_type;
	typedef GenotypeEncoderStatistics    stats_type;

public:
	GenotypeEncoder();
	GenotypeEncoder(const uint64_t samples);
	~GenotypeEncoder();

	inline void SetSamples(const uint64_t samples){ this->n_samples = samples; }
	inline const stats_type& GetUsageStats(void) const{ return(this->stats_); }

	bool Encode(const containers::VcfContainer& container, yon1_vnt_t* rcds, block_type& block, const yon_gt_ppa& ppa) const;
	bool Encode(yon1_vnt_t* rcds, const uint32_t n_rcds, block_type& block, const yon_gt_ppa& ppa) const;

	yon_gt_assess Assess(const bcf1_t* entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess Assess(const yon1_vnt_t& entry, const GenotypeSummary& gr_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessDiploidBiallelic(const bcf1_t* entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessDiploidBiallelic(const yon1_vnt_t& entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessDiploidMultiAllelic(const bcf1_t* entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessDiploidMultiAllelic(const yon1_vnt_t& entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessMultiploid(const bcf1_t* entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;
	yon_gt_assess AssessMultiploid(const yon1_vnt_t& entry, const GenotypeSummary& gt_summary, const yon_gt_ppa& ppa) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidBiallelic(const bcf1_t* entry,
                                    const GenotypeSummary& gt_summary,
	                                const yon_gt_ppa& ppa,
                                    container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidBiallelic(const yon1_vnt_t& entry,
									const GenotypeSummary& gt_summary,
									const yon_gt_ppa& ppa,
									container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidMultiAllelic(const bcf1_t* entry,
	                                   const GenotypeSummary& gt_summary,
	                                   const yon_gt_ppa& ppa,
	                                   container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidMultiAllelic(const yon1_vnt_t& entry,
									   const GenotypeSummary& gt_summary,
									   const yon_gt_ppa& ppa,
									   container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeMultiploid(const bcf1_t* entry,
	                          const GenotypeSummary& gt_summary,
	                          const yon_gt_ppa& ppa,
	                          container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeMultiploid(const yon1_vnt_t& entry,
							  const GenotypeSummary& gt_summary,
							  const yon_gt_ppa& ppa,
							  container_type& dst) const;

private:
	uint64_t n_samples; // number of samples
	stats_type stats_;
};

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeDiploidBiallelic(const bcf1_t* entry,
                                                 const GenotypeSummary& gt_summary,
                                                 const yon_gt_ppa& ppa,
                                                 container_type& dst) const
{
	assert(entry->d.fmt[0].n == 2);
	assert(entry->n_allele == 2);
	assert(gt_summary.n_vector_end == 0);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(ppa.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	uint64_t n_runs = 0; // Number of runs.
	uint64_t l_runs = 1; // Current run length.

	// 1 + hasMissing + hasMixedPhasing
	const uint8_t shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const uint8_t add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	const uint64_t limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (base_ploidy*shift + add)) - 1;
	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID(gt[ppa[0] * sizeof(int8_t) * base_ploidy],
	                                                   gt[ppa[0] * sizeof(int8_t) * base_ploidy + 1],
	                                                   shift, add);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[ppa[i] * sizeof(int8_t) * base_ploidy];
		uint32_t rle_ppa_current = YON_PACK_GT_DIPLOID(gt_ppa_target[0], gt_ppa_target[1], shift, add);

		if(rle_ppa_current != rle_ppa_current_ref || l_runs == limit){
			YON_RLE_TYPE RLE = l_runs;
			RLE <<= (base_ploidy*shift + add);
			RLE |= rle_ppa_current_ref;
			assert((RLE >> (base_ploidy*shift + add)) == l_runs);

			// Push RLE to buffer
			dst.AddLiteral((YON_RLE_TYPE)RLE);
			++dst.header.n_additions;

			rle_ppa_current_ref = rle_ppa_current;

			l_runs = 0;
			++n_runs;
		}
		++l_runs;
	}
	++n_runs;

	YON_RLE_TYPE RLE = l_runs;
	RLE <<= (base_ploidy*shift + add);
	RLE |= rle_ppa_current_ref;
	assert((RLE >> (base_ploidy*shift + add)) == l_runs);

	// Push RLE to buffer
	dst.AddLiteral((YON_RLE_TYPE)RLE);
	++dst.header.n_additions;
	++dst.header.n_entries;

	//std::cerr << n_runs << "\t" << dst.buffer_data_uncompressed.size() << std::endl;

	return n_runs;
}

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeDiploidBiallelic(const yon1_vnt_t& entry,
                                                 const GenotypeSummary& gt_summary,
                                                 const yon_gt_ppa& ppa,
                                                 container_type& dst) const
{
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	assert(entry.gt->m == 2);
	assert(entry.n_alleles == 2);
	assert(ppa.n_s == entry.gt->n_s);
	const uint8_t base_ploidy = entry.gt->m;

	uint64_t n_runs = 0; // Number of runs.
	uint64_t l_runs = 1; // Current run length.

	// 1 + hasMissing + hasMixedPhasing
	const uint8_t shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const uint8_t add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	const uint64_t limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (base_ploidy*shift + add)) - 1;
	uint32_t rle_ppa_current_ref = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[ppa[0]], shift, add);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		uint32_t rle_ppa_current = YON_PACK_GT_RCD_DIPLOID_EXPAND(entry.gt->d_exp[ppa[i]], shift, add);

		if(rle_ppa_current != rle_ppa_current_ref || l_runs == limit){
			YON_RLE_TYPE RLE = l_runs;
			RLE <<= (base_ploidy*shift + add);
			RLE |= rle_ppa_current_ref;
			assert((RLE >> (base_ploidy*shift + add)) == l_runs);

			// Push RLE to buffer
			dst.AddLiteral((YON_RLE_TYPE)RLE);
			++dst.header.n_additions;

			rle_ppa_current_ref = rle_ppa_current;

			l_runs = 0;
			++n_runs;
		}
		++l_runs;
	}
	++n_runs;

	YON_RLE_TYPE RLE = l_runs;
	RLE <<= (base_ploidy*shift + add);
	RLE |= rle_ppa_current_ref;
	assert((RLE >> (base_ploidy*shift + add)) == l_runs);

	// Push RLE to buffer
	dst.AddLiteral((YON_RLE_TYPE)RLE);
	++dst.header.n_additions;
	++dst.header.n_entries;

	//std::cerr << n_runs << "\t" << dst.buffer_data_uncompressed.size() << std::endl;

	return n_runs;
}

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeDiploidMultiAllelic(const bcf1_t* entry,
                                                    const GenotypeSummary& gt_summary,
                                                    const yon_gt_ppa& ppa,
                                                    container_type& dst) const
{
	assert(entry->d.fmt[0].n == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(ppa.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	int64_t n_runs = 0; // Number of runs.
	int64_t l_runs = 1; // Current run length.

	// Assess RLE cost
	const uint8_t shift = ceil(log2(entry->n_allele + 2 + 1));
	const uint8_t add   = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size.
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
	int64_t limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (base_ploidy*shift + add)) - 1;
	assert(limit > 0);

	uint8_t gt_remap[256];
	memset(gt_remap, 256, 255);
	for(uint32_t i = 0; i <= entry->n_allele; ++i){
		gt_remap[i << 1]       = ((i+1) << 1);
		gt_remap[(i << 1) + 1] = ((i+1) << 1) + 1;
	}
	gt_remap[0]   = 0;
	gt_remap[129] = 1;

	// Initial reference entry.
	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[ppa[0] * sizeof(int8_t) * base_ploidy]] >> 1,
														   gt_remap[gt[ppa[0] * sizeof(int8_t) * base_ploidy + 1]] >> 1,
														   shift, add,
														   gt_remap[gt[ppa[0] * sizeof(int8_t) * base_ploidy + 1]]);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[ppa[i] * sizeof(int8_t) * base_ploidy];
		uint32_t rle_ppa_current = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt_ppa_target[0]] >> 1,
														   gt_remap[gt_ppa_target[1]] >> 1,
														   shift, add,
														   gt_remap[gt_ppa_target[1]]);

		assert(gt_remap[gt_ppa_target[0]] != 255);
		assert(gt_remap[gt_ppa_target[1]] != 255);

		if(rle_ppa_current != rle_ppa_current_ref || l_runs == limit){
			YON_RLE_TYPE RLE = l_runs;
			RLE <<= (base_ploidy*shift + add);
			RLE |= rle_ppa_current_ref;
			assert((RLE >> (base_ploidy*shift + add)) == l_runs);

			// Push RLE to buffer
			dst.AddLiteral((YON_RLE_TYPE)RLE);
			++dst.header.n_additions;

			rle_ppa_current_ref = rle_ppa_current;

			l_runs = 0;
			++n_runs;
		}

		++l_runs;
	}
	++n_runs;

	YON_RLE_TYPE RLE = l_runs;
	RLE <<= (base_ploidy*shift + add);
	RLE |= rle_ppa_current_ref;
	assert((RLE >> (base_ploidy*shift + add)) == l_runs);

	// Push RLE to buffer
	dst.AddLiteral((YON_RLE_TYPE)RLE);
	++dst.header.n_additions;
	++dst.header.n_entries;

	return n_runs;

	//std::cerr << n_runs << "\t" << dst.buffer_data_uncompressed.size() << std::endl;
}

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeDiploidMultiAllelic(const yon1_vnt_t& entry,
                                                    const GenotypeSummary& gt_summary,
                                                    const yon_gt_ppa& ppa,
                                                    container_type& dst) const
{
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	assert(entry.gt->m == 2);
	assert(ppa.n_s == entry.gt->n_s);
	const uint8_t base_ploidy = entry.gt->m;

	int64_t n_runs = 0; // Number of runs.
	int64_t l_runs = 1; // Current run length.

	// Assess RLE cost
	const uint8_t shift = ceil(log2(entry.n_alleles + 2 + 1));
	const uint8_t add   = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	// Values set to signed integers as values can underflow if
	// the do not fit in the word size.
	// Ploidy*shift_size bits for alleles and 1 bit for phase information (if required)
	// Cost: 2^(8*word_width - (ploidy*(n_alleles + has_missing + hasEOV + 1) + has_mixed_phasing))
	int64_t limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (base_ploidy*shift + add)) - 1;
	assert(limit > 0);

	// Initial reference entry.
	uint32_t rle_ppa_current_ref = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[ppa[0]], shift, add);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		uint32_t rle_ppa_current = YON_PACK_GT_RCD_NALLELIC_EXPAND(entry.gt->d_exp[ppa[i]], shift, add);

		if(rle_ppa_current != rle_ppa_current_ref || l_runs == limit){
			YON_RLE_TYPE RLE = l_runs;
			RLE <<= (base_ploidy*shift + add);
			RLE |= rle_ppa_current_ref;
			assert((RLE >> (base_ploidy*shift + add)) == l_runs);

			// Push RLE to buffer
			dst.AddLiteral((YON_RLE_TYPE)RLE);
			++dst.header.n_additions;

			rle_ppa_current_ref = rle_ppa_current;

			l_runs = 0;
			++n_runs;
		}

		++l_runs;
	}
	++n_runs;

	YON_RLE_TYPE RLE = l_runs;
	RLE <<= (base_ploidy*shift + add);
	RLE |= rle_ppa_current_ref;
	assert((RLE >> (base_ploidy*shift + add)) == l_runs);

	// Push RLE to buffer
	dst.AddLiteral((YON_RLE_TYPE)RLE);
	++dst.header.n_additions;
	++dst.header.n_entries;

	return n_runs;

	//std::cerr << n_runs << "\t" << dst.buffer_data_uncompressed.size() << std::endl;
}

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeMultiploid(const bcf1_t* entry,
                                           const GenotypeSummary& gt_summary,
                                           const yon_gt_ppa& ppa,
                                           container_type& dst) const
{
	// This method is currently only valid if the genotypic
	// data is stored as BCF_BT_INT8 in the htslib bcf1_t
	// record.
	assert(entry->d.fmt[0].type == BCF_BT_INT8);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	const uint64_t limit = std::numeric_limits<YON_RLE_TYPE>::max();
	assert(ppa.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	// Remap genotype encoding such that 0 maps to missing and
	// 1 maps to the sentinel node symbol (EOV).
	uint8_t gt_remap[256];
	memset(gt_remap, 256, 255);
	for(uint32_t i = 0; i <= entry->n_allele; ++i){
		gt_remap[i << 1]       = ((i+1) << 1);
		gt_remap[(i << 1) + 1] = ((i+1) << 1) + 1;
	}
	gt_remap[0]   = 0;
	gt_remap[129] = 1;

	// Start parameters for run-length encoding.
	uint64_t n_runs = 0;     // Number of runs.
	YON_RLE_TYPE l_run  = 1; // Current run length.

	// Keep track of the current reference sequence as we
	// iterate over the available genotypes.
	uint8_t* reference = new uint8_t[base_ploidy];
	for(uint32_t i = 0; i < base_ploidy; ++i)
		reference[i] = gt_remap[gt[ppa[0] * sizeof(int8_t) * base_ploidy + i]];

	// Hash of current reference genotype sequence.
	uint64_t hash_value_ppa_ref = XXH64(&gt[ppa[0] * sizeof(int8_t) * base_ploidy], sizeof(int8_t) * base_ploidy, 89231478);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[ppa[i] * sizeof(int8_t) * base_ploidy];
		uint64_t hash_value_ppa = XXH64(gt_ppa_target, sizeof(int8_t) * base_ploidy, 89231478);

		if(hash_value_ppa != hash_value_ppa_ref){
			dst.AddLiteral(l_run);
			for(uint32_t k = 0; k < base_ploidy; ++k) dst.AddLiteral(reference[k]);
			++dst.header.n_additions;

			++n_runs;
			l_run = 0;
			hash_value_ppa_ref = hash_value_ppa;
			for(uint32_t k = 0; k < base_ploidy; ++k) reference[k] = gt_remap[gt_ppa_target[k]];
		}

		// Overflow: trigger a break
		if(l_run == limit){
			dst.AddLiteral(l_run);
			for(uint32_t k = 0; k < base_ploidy; ++k) dst.AddLiteral(reference[k]);
			++dst.header.n_additions;

			++n_runs;
			l_run = 0;
		}
		++l_run;
	}
	dst.AddLiteral(l_run);
	for(uint32_t k = 0; k < base_ploidy; ++k) dst.AddLiteral(reference[k]);
	++dst.header.n_additions;
	++dst.header.n_entries;

	//std::cerr << l_run << ":" << (uint32_t)reference[0];
	//for(uint32_t k = 1; k < base_ploidy; ++k) std::cerr << "," << (uint32_t)reference[k];
	//std::cerr << std::endl;
	++n_runs;

	//std::cerr << n_runs << "\t" << dst.buffer_data_uncompressed.size() << std::endl;

	delete [] reference;

	return n_runs;
}

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeMultiploid(const yon1_vnt_t& entry,
                                           const GenotypeSummary& gt_summary,
                                           const yon_gt_ppa& ppa,
                                           container_type& dst) const
{
	// This method is currently only valid if the genotypic
	// data is stored as BCF_BT_INT8 in the htslib bcf1_t
	// record.
	assert(entry.gt != nullptr);
	assert(entry.gt->d_exp != nullptr);
	const uint8_t   base_ploidy = entry.gt->m;
	assert(ppa.n_s == entry.gt->n_s);

	const uint64_t limit = std::numeric_limits<YON_RLE_TYPE>::max();

	std::cerr << "not implemented yet" << std::endl;
	exit(1);

	return 0;
}

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
