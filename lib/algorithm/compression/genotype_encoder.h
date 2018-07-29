#ifndef ENCODERGENOTYPESRLE_H_
#define ENCODERGENOTYPESRLE_H_

#include <algorithm>
#include <bitset>
#include <cassert>
#include <thread>

#include "core/genotype_summary.h"
#include "io/bcf/bcf_reader.h"
#include "containers/variant_block.h"
#include "core/variant_controller.h"

#include "io/htslib_integration.h"
#include "algorithm/permutation/radix_sort_gt.h"

namespace tachyon{
namespace algorithm{

#define ENCODER_GT_DEBUG 0
#define YON_PACK_GT_DIPLOID(A, B, SHIFT, ADD)                 (bcf::BCF_UNPACK_GENOTYPE(A) << ((SHIFT) + (ADD))) | (bcf::BCF_UNPACK_GENOTYPE(B) << (ADD)) | ((A) & (ADD))
#define YON_PACK_GT_DIPLOID_NALLELIC(A, B, SHIFT, ADD, PHASE) ((A) << ((SHIFT) + (ADD))) | ((B) << (ADD)) | ((PHASE) & (ADD))

struct GenotypeEncoderStatistics{
	GenotypeEncoderStatistics(){
		memset(this->rle_counts,         0, sizeof(U64)*4);
		memset(this->rle_simple_counts,  0, sizeof(U64)*4);
		memset(this->diploid_bcf_counts, 0, sizeof(U64)*3);
		memset(this->bcf_counts,         0, sizeof(U64)*3);
	}

	U64 getTotal(void) const{
		U64 total = 0;
		for(U32 i = 0; i < 4; ++i) total += this->rle_counts[i];
		for(U32 i = 0; i < 4; ++i) total += this->rle_simple_counts[i];
		for(U32 i = 0; i < 3; ++i) total += this->diploid_bcf_counts[i];
		for(U32 i = 0; i < 3; ++i) total += this->bcf_counts[i];
		return(total);
	}

	U64 rle_counts[4];
	U64 rle_simple_counts[4];
	U64 diploid_bcf_counts[3];
	U64 bcf_counts[3];
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

	GenotypeEncoderSlaveHelper(const U32 start_capacity) :
		encoding_type(YON_GT_RLE_DIPLOID_BIALLELIC),
		gt_primitive(YON_GT_BYTE),
		n_runs(0)
	{
		// We only use the uncompressed buffer
		// no strides or compressed buffers
		container.buffer_data_uncompressed.resize(start_capacity);
	}
	~GenotypeEncoderSlaveHelper(){}

	// Overload operator += for block and RTYPE helper
	friend block_type& operator+=(block_type& block, const self_type& helper){
		block.base_containers[YON_BLK_GT_SUPPORT].Add((U32)helper.n_runs);
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
	U32 n_runs;
	container_type container;
};

class GenotypeEncoder {
private:
	typedef GenotypeEncoder              self_type;
	typedef io::BasicBuffer              buffer_type;
	typedef bcf::BCFReader               bcf_reader_type;
	typedef bcf::BCFEntry                bcf_type;
	typedef core::MetaEntry              meta_type;
	typedef containers::DataContainer    container_type;
	typedef containers::VariantBlock     block_type;
	typedef GenotypeEncoderStatistics    stats_type;

	typedef struct __RLEAssessHelper{
		explicit __RLEAssessHelper(void) :
				word_width(1),
				n_runs(0)
		{}
		__RLEAssessHelper(const BYTE& word_width,
				          const U64& n_runs) :
			word_width(word_width),
			n_runs(n_runs)
		{}
		~__RLEAssessHelper(){}

		BYTE word_width;
		U64 n_runs;

	} rle_helper_type;

public:
	GenotypeEncoder();
	GenotypeEncoder(const U64 samples);
	~GenotypeEncoder();
	bool Encode(const bcf_type& bcf_entry, meta_type& meta, block_type& block, const U32* const ppa);
	bool EncodeParallel(const bcf_reader_type& bcf_reader, meta_type* meta_entries, block_type& block, const U32* const ppa, const U32 n_threads);
	bool EncodeParallel(const bcf_type& bcf_entry, meta_type& meta, const U32* const ppa, GenotypeEncoderSlaveHelper& slave_helper) const;
	inline void setSamples(const U64 samples){ this->n_samples = samples; }
	inline const stats_type& getUsageStats(void) const{ return(this->stats_); }

	bool Assess(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const;
	bool AssessHaploid(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const;
	bool AssessDiploidBiallelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const;
	bool AssessDiploidMultiAllelic(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const;
	bool AssessMultiploid(const bcf1_t* entry, const io::VcfGenotypeSummary& gt_summary, const algorithm::yon_gt_ppa& permutation_array) const;

private:
	const rle_helper_type assessDiploidRLEBiallelic(const bcf_type& bcf_entry, const U32* const ppa) const;
	const rle_helper_type assessDiploidRLEnAllelic(const bcf_type& bcf_entry, const U32* const ppa) const;
	//const rle_helper_type assessMploidRLEBiallelic(const bcf_type& bcf_entry, const U32* const ppa) const;
	//const rle_helper_type assessMploidRLEnAllelic(const bcf_type& bcf_entry, const U32* const ppa) const;

	template <class YON_STORE_TYPE, class BCF_GT_TYPE = BYTE> bool EncodeBCFStyle(const bcf_type& bcf_entry, container_type& container, U64& n_runs) const;
	template <class YON_RLE_TYPE, class BCF_GT_TYPE = BYTE> bool EncodeDiploidBCF(const bcf_type& bcf_entry, container_type& runs, U64& n_runs, const U32* const ppa) const;
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEBiallelic(const bcf_type& bcf_entry, container_type& runs, const U32* const ppa, const rle_helper_type& helper) const;
	template <class YON_RLE_TYPE> bool EncodeDiploidRLEnAllelic(const bcf_type& bcf_entry, container_type& runs, const U32* const ppa, const rle_helper_type& helper) const;
	//template <class T> bool EncodeMploidRLEBiallelic(const bcf_type& bcf_entry, container_type& runs, U64& n_runs, const U32* const ppa) const;
	//template <class T> bool EncodeMploidRLENallelic(const bcf_type& bcf_entry, container_type& runs, U64& n_runs, const U32* const ppa) const;

	/**<
	 * Supportive reduce function for updating local import statistics
	 * following parallel execution of `EncodeParallel`. Iteratively
	 * call this function with all subproblems to calculate the total
	 * import statistics of genotypes.
	 * @param helper Input helper structure
	 */
	void updateStatistics(const GenotypeEncoderSlaveHelper& helper);

private:
	U64 n_samples; // number of samples
	stats_type stats_;
};

template <class YON_STORE_TYPE, class BCF_GT_TYPE>
bool GenotypeEncoder::EncodeBCFStyle(const bcf_type& bcf_entry,
                                     container_type& simple,
                                                U64& n_runs) const
{
	const BYTE ploidy = bcf_entry.gt_support.ploidy;
	U32 bcf_gt_pos = bcf_entry.formatID[0].l_offset;
	const BCF_GT_TYPE missing_value = (BCF_GT_TYPE)1 << (sizeof(BCF_GT_TYPE)*8 - 1);
	const BCF_GT_TYPE EOV_value     = missing_value + 1;

	// Pack genotypes as
	// allele | phasing
	U32 j = 0;
	for(U32 i = 0; i < this->n_samples * ploidy; i += ploidy, ++j){
		for(U32 p = 0; p < ploidy; ++p){
			const BCF_GT_TYPE& allele = *reinterpret_cast<const BCF_GT_TYPE* const>(&bcf_entry.data[bcf_gt_pos]);
			if((allele >> 1) == 0) simple.AddLiteral((YON_STORE_TYPE)0); // missing
			else if(allele == EOV_value) simple.AddLiteral((YON_STORE_TYPE)1); // eov
			else { // otherwise
				// Add 1 because 1 is reserved for EOV
				const YON_STORE_TYPE val = ((allele >> 1) + 1) << 1 | (allele & 1);
				simple.AddLiteral((YON_STORE_TYPE)val);
			}
			bcf_gt_pos += sizeof(BCF_GT_TYPE);
		}
	}

	n_runs = this->n_samples*ploidy;
	simple.header.n_additions += n_runs;

	return(true);
}

template <class YON_RLE_TYPE, class BCF_GT_TYPE>
bool GenotypeEncoder::EncodeDiploidBCF(const bcf_type& bcf_entry,
		                               container_type& simple,
										          U64& n_runs,
                                     const U32* const  ppa) const
{
	const BYTE ploidy = 2;
	// Shift size is equivalent to floor((sizeof(T)*8 - 1)/2)
	const BYTE shift_size = (sizeof(YON_RLE_TYPE)*8 - 1) / 2;

	// Start of GT byte stream
	const char* const data = &bcf_entry.data[bcf_entry.formatID[0].l_offset];

	const BCF_GT_TYPE missing_value = (BCF_GT_TYPE)1 << (sizeof(BCF_GT_TYPE)*8 - 1);
	const BCF_GT_TYPE EOV_value     = missing_value + 1;
	// Pack genotypes as
	// allele A | allele B | phasing information
	U32 ppa_pos = 0;
	//YON_RLE_TYPE temp = 0;
	for(U32 i = 0; i < this->n_samples * ploidy; i += ploidy){
		BCF_GT_TYPE allele1 = *reinterpret_cast<const BCF_GT_TYPE* const>(&data[ploidy*sizeof(BCF_GT_TYPE)*ppa[ppa_pos]]);
		BCF_GT_TYPE allele2 = *reinterpret_cast<const BCF_GT_TYPE* const>(&data[ploidy*sizeof(BCF_GT_TYPE)*ppa[ppa_pos] + sizeof(BCF_GT_TYPE)]);
		const bool phasing = allele2 & 1;

		if((allele1 >> 1) == 0)       allele1 = 0;
		else if(allele1 == EOV_value) allele1 = 1;
		else allele1 = (allele1 >> 1) + 1;
		if((allele2 >> 1) == 0)       allele2 = 0;
		else if(allele2 == EOV_value) allele2 = 1;
		else allele2 = (allele2 >> 1) + 1;

		const YON_RLE_TYPE packed = (allele1 << (shift_size + 1)) |
				                    (allele2 << 1) |
									(phasing &  1);

		simple.AddLiteral((YON_RLE_TYPE)packed);
		++ppa_pos;
	}

	n_runs = this->n_samples;
	simple.header.n_additions += n_runs;

	return(true);
}

template <class YON_RLE_TYPE>
bool GenotypeEncoder::EncodeDiploidRLEBiallelic(const bcf_type& bcf_entry,
		                                        container_type& runs,
                                              const U32* const  ppa,
								         const rle_helper_type& helper) const
{
	const BYTE ploidy   = 2;
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	const BYTE shift    = bcf_entry.gt_support.hasMissing    ? 2 : 1;
	const BYTE add      = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;

	// Run limits
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE) - (ploidy*shift + add)) - 1;
	const char* const data       = &bcf_entry.data[bcf_entry.formatID[0].l_offset];
	const BYTE& allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0]]);
	const BYTE& allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[0] + sizeof(BYTE)]);
	YON_RLE_TYPE packed = YON_PACK_GT_DIPLOID(allele2, allele1, shift, add);

	U32 ppa_pos = 1;
	U64 n_runs = 0;
	for(U32 i = ploidy; i < this->n_samples * ploidy; i += ploidy){
		const BYTE& allele1 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos]]);
		const BYTE& allele2 = *reinterpret_cast<const BYTE* const>(&data[ploidy*sizeof(BYTE)*ppa[ppa_pos] + sizeof(BYTE)]);
		const YON_RLE_TYPE packed_internal = YON_PACK_GT_DIPLOID(allele2, allele1, shift, add);

		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = length;
			RLE <<= (ploidy*shift + add);
			RLE |= packed;
			assert((RLE >> (ploidy*shift + add)) == length);

			// Push RLE to buffer
			runs.AddLiteral((YON_RLE_TYPE)RLE);

			// Reset and update
			sumLength += length;
			length = 0;
			packed = packed_internal;
			++n_runs;
		}
		++length;
		++ppa_pos;
	}
	// Last entry
	// Prepare RLE
	RLE = length;
	RLE <<= (ploidy*shift + add);
	RLE |= packed;
	assert((RLE >> (ploidy*shift + add)) == length);

	// Push RLE to buffer
	runs.AddLiteral((YON_RLE_TYPE)RLE);
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	assert(helper.n_runs == n_runs);
	runs.header.n_additions += n_runs;
	assert(ppa_pos == n_samples);

#if ENCODER_GT_DEBUG == 1
	std::cout << 0 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif
	return(true);
}

template <class YON_RLE_TYPE>
bool GenotypeEncoder::EncodeDiploidRLEnAllelic(const bcf_type& bcf_entry,
		                                       container_type& runs,
											 const U32* const  ppa,
										const rle_helper_type& helper) const
{
	const BYTE ploidy   = 2;
	U32 sumLength       = 0;
	YON_RLE_TYPE length = 1;
	YON_RLE_TYPE RLE    = 0;
	const BYTE shift    = ceil(log2(bcf_entry.body->n_allele + 2 + 1)); // Bits occupied per allele
	const BYTE add      = bcf_entry.gt_support.mixedPhasing  ? 1 : 0;
	const YON_RLE_TYPE run_limit = pow(2, 8*sizeof(YON_RLE_TYPE)  - (ploidy*shift + add)) - 1;

	// Setup first run
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

	YON_RLE_TYPE packed  = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add, phase);

	U32 ppa_pos = 1;
	U64 n_runs = 0;
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

		const YON_RLE_TYPE packed_internal = YON_PACK_GT_DIPLOID_NALLELIC(allele2, allele1, shift, add, phase);


		if(packed != packed_internal || length == run_limit){
			// Prepare RLE
			RLE = (length << (ploidy*shift + add)) | packed;
			assert((RLE >> (ploidy*shift + add)) == length);
			assert(length != 0);

			// Push RLE to buffer
			runs.AddLiteral((YON_RLE_TYPE)RLE);

			// Reset and update
			sumLength += length;
			length = 0;
			packed = packed_internal;
			++n_runs;
		}
		++length;
		++ppa_pos;
	}
	// Last entry
	// Prepare RLE
	RLE = (length << (ploidy*shift + add)) | packed;
	assert((RLE >> (ploidy*shift + add)) == length);
	assert(length != 0);

	// Push RLE to buffer
	runs.AddLiteral((YON_RLE_TYPE)RLE);
	++n_runs;

	// Reset and update
	sumLength += length;
	assert(sumLength == this->n_samples);
	assert(helper.n_runs == n_runs);
	assert(ppa_pos == n_samples);
	runs.header.n_additions += n_runs;

#if ENCODER_GT_DEBUG == 1
	std::cout << 1 << '\t' << n_runs << '\t' << sizeof(YON_RLE_TYPE) << '\n';
#endif

	return(true);
}

/**<
 * Parallel support structure: this object encapsulates
 * a thread that runs the `EncodeParallel` function with
 * a stride size of N_THREADS
 */
struct CalcSlave{
	typedef CalcSlave       self_type;
	typedef bcf::BCFEntry   bcf_type;
	typedef bcf::BCFReader  bcf_reader_type;
	typedef core::MetaEntry meta_type;
	typedef GenotypeEncoderSlaveHelper helper_type;

	CalcSlave() :
		thread_idx_(0),
		n_threads_(0),
		encoder_(nullptr),
		reader_(nullptr),
		meta_entries_(nullptr),
		ppa_(nullptr),
		helpers_(nullptr)
	{}

	~CalcSlave(){}

	std::thread* Start(const GenotypeEncoder& encoder,
		               const U32 thread_idx,
		               const U32 n_threads,
		               const bcf_reader_type& reader,
		               meta_type* meta_entries,
		               const U32* const ppa,
		               helper_type* helpers)
	{
		this->encoder_      = &encoder;
		this->thread_idx_   = thread_idx;
		this->n_threads_    = n_threads;
		this->reader_       = &reader;
		this->meta_entries_ = meta_entries;
		this->ppa_          = ppa;
		this->helpers_      = helpers;

		this->thread = std::thread(&self_type::Run_, this);
		return(&this->thread);
	}

private:
	void Run_(void){
		for(U32 i = this->thread_idx_; i < this->reader_->size(); i += this->n_threads_){
			this->encoder_->EncodeParallel((*this->reader_)[i], this->meta_entries_[i], this->ppa_, this->helpers_[i]);
		}
	}

private:
	U32 thread_idx_;
	U32 n_threads_;
	const GenotypeEncoder* encoder_;
	const bcf_reader_type* reader_;
	meta_type* meta_entries_;
	const U32* ppa_;
	helper_type* helpers_;

public:
	std::thread thread;
};

}
}

#endif /* ENCODERGENOTYPESRLE_H_ */
