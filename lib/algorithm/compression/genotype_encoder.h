#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>
#include <thread>

#include "containers/variant_block.h"
#include "core/variant_controller.h"
#include "core/genotypes.h"
#include "io/vcf_utils.h"
#include "containers/vcf_container.h"

namespace tachyon{
namespace algorithm{

const uint8_t BCF_UNPACK_TACHYON[3] = {2, 0, 1};
#define BCF_UNPACK_GENOTYPE(A) BCF_UNPACK_TACHYON[((A) >> 1)]
const char BCF_TYPE_SIZE[8] = {0,1,2,4,0,4,0,1};

#define ENCODER_GT_DEBUG 0
#define YON_PACK_GT_DIPLOID(A, B, SHIFT, ADD)                 (BCF_UNPACK_GENOTYPE(A) << ((SHIFT) + (ADD))) | (BCF_UNPACK_GENOTYPE(B) << (ADD)) | ((A) & (ADD))
#define YON_PACK_GT_DIPLOID_NALLELIC(A, B, SHIFT, ADD, PHASE) ((A) << ((SHIFT) + (ADD))) | ((B) << (ADD)) | ((PHASE) & (ADD))

struct yon_gt_assess {
	uint8_t GetCheapestPrimitive(void) const{
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

/**<
 * Supportive structure for the parallel import of genotypes
 */
struct GenotypeEncoderSlaveHelper{
	typedef GenotypeEncoderSlaveHelper self_type;
	typedef containers::DataContainer  container_type;
	typedef containers::VariantBlock   block_type;

public:
	GenotypeEncoderSlaveHelper() :
		encoding_type(YON_GT_RLE_DIPLOID_BIALLELIC),
		gt_primitive(YON_GT_BYTE),
		n_runs(0)
	{
	}

	GenotypeEncoderSlaveHelper(const uint32_t start_capacity) :
		encoding_type(YON_GT_RLE_DIPLOID_BIALLELIC),
		gt_primitive(YON_GT_BYTE),
		n_runs(0)
	{
		// We only use the uncompressed buffer
		// no strides or compressed buffers
		container.data_uncompressed.resize(start_capacity);
	}
	~GenotypeEncoderSlaveHelper(){}

	// Overload operator += for block and RTYPE helper
	friend block_type& operator+=(block_type& block, const self_type& helper){
		block.base_containers[YON_BLK_GT_SUPPORT].Add((uint32_t)helper.n_runs);
		++block.base_containers[YON_BLK_GT_SUPPORT];

		if(helper.encoding_type == YON_GT_RLE_DIPLOID_BIALLELIC){
			if(helper.gt_primitive == YON_GT_BYTE){
				block.base_containers[YON_BLK_GT_INT8] += helper.container;
				++block.base_containers[YON_BLK_GT_INT8];
			} else if(helper.gt_primitive == YON_GT_U16){
				block.base_containers[YON_BLK_GT_INT16] += helper.container;
				++block.base_containers[YON_BLK_GT_INT16];
			} else if(helper.gt_primitive == YON_GT_U32){
				block.base_containers[YON_BLK_GT_INT32] += helper.container;
				++block.base_containers[YON_BLK_GT_INT32];
			} else if(helper.gt_primitive == YON_GT_U64){
				block.base_containers[YON_BLK_GT_INT64] += helper.container;
				++block.base_containers[YON_BLK_GT_INT64];
			}
		} else if(helper.encoding_type == YON_GT_RLE_DIPLOID_NALLELIC){
			if(helper.gt_primitive == YON_GT_BYTE){
				block.base_containers[YON_BLK_GT_S_INT8] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT8];
			} else if(helper.gt_primitive == YON_GT_U16){
				block.base_containers[YON_BLK_GT_S_INT16] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT16];
			} else if(helper.gt_primitive == YON_GT_U32){
				block.base_containers[YON_BLK_GT_S_INT32] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT32];
			} else if(helper.gt_primitive == YON_GT_U64){
				block.base_containers[YON_BLK_GT_S_INT64] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT64];
			}
		} else if(helper.encoding_type == YON_GT_BCF_DIPLOID){
			if(helper.gt_primitive == YON_GT_BYTE){
				block.base_containers[YON_BLK_GT_S_INT8] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT8];
			} else if(helper.gt_primitive == YON_GT_U16){
				block.base_containers[YON_BLK_GT_S_INT16] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT16];
			} else if(helper.gt_primitive == YON_GT_U32){
				block.base_containers[YON_BLK_GT_S_INT32] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT32];
			}
		} else if(helper.encoding_type == YON_GT_BCF_STYLE){
			if(helper.gt_primitive == YON_GT_BYTE){
				block.base_containers[YON_BLK_GT_S_INT8] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT8];
			} else if(helper.gt_primitive == YON_GT_U16){
				block.base_containers[YON_BLK_GT_S_INT16] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT16];
			} else if(helper.gt_primitive == YON_GT_U32){
				block.base_containers[YON_BLK_GT_S_INT32] += helper.container;
				++block.base_containers[YON_BLK_GT_S_INT32];
			}
		}

		return(block);
	}

public:
	TACHYON_GT_ENCODING       encoding_type;
	TACHYON_GT_PRIMITIVE_TYPE gt_primitive;
	uint32_t n_runs;
	container_type container;
};

class GenotypeEncoder {
public:
	typedef GenotypeEncoder              self_type;
	typedef io::BasicBuffer              buffer_type;
	typedef core::MetaEntry              meta_type;
	typedef containers::DataContainer    container_type;
	typedef containers::VariantBlock     block_type;
	typedef GenotypeEncoderStatistics    stats_type;

	typedef struct __RLEAssessHelper{
		explicit __RLEAssessHelper(void) :
				word_width(1),
				n_runs(0)
		{}
		__RLEAssessHelper(const uint8_t& word_width,
				          const uint64_t& n_runs) :
			word_width(word_width),
			n_runs(n_runs)
		{}
		~__RLEAssessHelper(){}

		uint8_t word_width;
		uint64_t n_runs;

	} rle_helper_type;

public:
	GenotypeEncoder();
	GenotypeEncoder(const uint64_t samples);
	~GenotypeEncoder();

	inline void SetSamples(const uint64_t samples){ this->n_samples = samples; }
	inline const stats_type& GetUsageStats(void) const{ return(this->stats_); }

	bool Encode(const containers::VcfContainer& container, meta_type* meta_entries, block_type& block, const yon_gt_ppa& permutation_array) const;
	bool EncodeParallel(const containers::VcfContainer& container, meta_type* meta_entries, GenotypeEncoderSlaveHelper& slave_helper, const yon_gt_ppa& permutation_array) const;

	yon_gt_assess Assess(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const yon_gt_ppa& permutation_array) const;
	yon_gt_assess AssessDiploidBiallelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const yon_gt_ppa& permutation_array) const;
	yon_gt_assess AssessDiploidMultiAllelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const yon_gt_ppa& permutation_array) const;
	yon_gt_assess AssessMultiploid(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const yon_gt_ppa& permutation_array) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidBiallelic(const bcf1_t* entry,
                                    const io::VcfGenotypeSummary& gt_summary,
	                                const yon_gt_ppa& permutation_array,
                                    container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeDiploidMultiAllelic(const bcf1_t* entry,
	                                   const io::VcfGenotypeSummary& gt_summary,
	                                   const yon_gt_ppa& permutation_array,
	                                   container_type& dst) const;

	template <class YON_RLE_TYPE>
	uint64_t EncodeMultiploid(const bcf1_t* entry,
	                          const io::VcfGenotypeSummary& gt_summary,
	                          const yon_gt_ppa& permutation_array,
	                          container_type& dst) const;

private:
	/**<
	 * Supportive reduce function for updating local import statistics
	 * following parallel execution of `EncodeParallel`. Iteratively
	 * call this function with all subproblems to calculate the total
	 * import statistics of genotypes.
	 * @param helper Input helper structure
	 */
	void updateStatistics(const GenotypeEncoderSlaveHelper& helper);

private:
	uint64_t n_samples; // number of samples
	stats_type stats_;
};

template <class YON_RLE_TYPE>
uint64_t GenotypeEncoder::EncodeDiploidBiallelic(const bcf1_t* entry,
                                                 const io::VcfGenotypeSummary& gt_summary,
                                                 const yon_gt_ppa& permutation_array,
                                                 container_type& dst) const
{
	assert(entry->d.fmt[0].n == 2);
	assert(entry->n_allele == 2);
	assert(gt_summary.n_vector_end == 0);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_s * base_ploidy == l_gt);
	assert(base_ploidy * sizeof(uint8_t) * this->n_samples == l_gt);

	uint64_t n_runs = 0; // Number of runs.
	uint64_t l_runs = 1; // Current run length.

	// 1 + hasMissing + hasMixedPhasing
	const uint8_t shift  = gt_summary.n_missing      ? 2 : 1; // 1-bits enough when no data missing {0,1}, 2-bits required when missing is available {0,1,2}
	const uint8_t add    = gt_summary.mixed_phasing  ? 1 : 0;

	// Run limits
	const uint64_t limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (base_ploidy*shift + add)) - 1;
	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID(gt[permutation_array[0] * sizeof(int8_t) * base_ploidy],
												  gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1],
												  shift, add);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
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
uint64_t GenotypeEncoder::EncodeDiploidMultiAllelic(const bcf1_t* entry,
                                                    const io::VcfGenotypeSummary& gt_summary,
                                                    const yon_gt_ppa& permutation_array,
                                                    container_type& dst) const
{
	assert(entry->d.fmt[0].n == 2);
	const uint8_t   base_ploidy = entry->d.fmt[0].n;
	const uint8_t*  gt   = entry->d.fmt[0].p;
	const uint32_t  l_gt = entry->d.fmt[0].p_len;
	assert(permutation_array.n_s * base_ploidy == l_gt);
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
	uint32_t rle_ppa_current_ref = YON_PACK_GT_DIPLOID_NALLELIC(gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy]] >> 1,
														   gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]] >> 1,
														   shift, add,
														   gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + 1]]);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
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
uint64_t GenotypeEncoder::EncodeMultiploid(const bcf1_t* entry,
                                           const io::VcfGenotypeSummary& gt_summary,
                                           const yon_gt_ppa& permutation_array,
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
	assert(permutation_array.n_s * base_ploidy == l_gt);
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
		reference[i] = gt_remap[gt[permutation_array[0] * sizeof(int8_t) * base_ploidy + i]];

	// Hash of current reference genotype sequence.
	uint64_t hash_value_ppa_ref = XXH64(&gt[permutation_array[0] * sizeof(int8_t) * base_ploidy], sizeof(int8_t) * base_ploidy, 89231478);

	// Iterate over all available samples.
	for(uint32_t i = 1; i < this->n_samples; ++i){
		const uint8_t* gt_ppa_target = &gt[permutation_array[i] * sizeof(int8_t) * base_ploidy];
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


/**<
 * Parallel support structure: this object encapsulates
 * a thread that runs the `EncodeParallel` function with
 * a stride size of N_THREADS
 */
struct CalcSlave{
	typedef CalcSlave       self_type;
	typedef core::MetaEntry meta_type;
	typedef GenotypeEncoderSlaveHelper helper_type;

	CalcSlave() :
		thread_idx_(0),
		n_threads_(0),
		encoder_(nullptr),
		//reader_(nullptr),
		meta_entries_(nullptr),
		ppa_(nullptr),
		helpers_(nullptr)
	{}

	~CalcSlave(){}

	std::thread* Start(const GenotypeEncoder& encoder,
		               const uint32_t thread_idx,
		               const uint32_t n_threads,
		               //const bcf_reader_type& reader,
		               meta_type* meta_entries,
		               const uint32_t* const ppa,
		               helper_type* helpers)
	{
		this->encoder_      = &encoder;
		this->thread_idx_   = thread_idx;
		this->n_threads_    = n_threads;
		//this->reader_       = &reader;
		this->meta_entries_ = meta_entries;
		this->ppa_          = ppa;
		this->helpers_      = helpers;

		this->thread = std::thread(&self_type::Run_, this);
		return(&this->thread);
	}

private:
	void Run_(void){
		exit(1);
		for(uint32_t i = this->thread_idx_; i < 5; i += this->n_threads_){
			//this->encoder_->EncodeParallel((*this->reader_)[i], this->meta_entries_[i], this->ppa_, this->helpers_[i]);
		}
	}

private:
	uint32_t thread_idx_;
	uint32_t n_threads_;
	const GenotypeEncoder* encoder_;
	//const bcf_reader_type* reader_;
	meta_type* meta_entries_;
	const uint32_t* ppa_;
	helper_type* helpers_;

public:
	std::thread thread;
};

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
