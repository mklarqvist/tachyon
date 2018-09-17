#ifndef CORE_TS_TV_OBJECT_H_
#define CORE_TS_TV_OBJECT_H_

#include "variant_record.h"

namespace tachyon{

#define YON_GT_TSTV_A       0
#define YON_GT_TSTV_T       1
#define YON_GT_TSTV_G       2
#define YON_GT_TSTV_C       3
#define YON_GT_TSTV_UNKNOWN 4
#define YON_GT_TSTV_MISS    5
#define YON_GT_TSTV_EOV     6
#define YON_GT_TSTV_INS     7
#define YON_GT_TSTV_DEL     8

// ASCII positions for A,T,G,C are filled with their respective
// encodings and the background is filled with YON_GT_TSTV_UNKNOWN.
const uint8_t YON_STATS_TSTV_LOOKUP[256] =
{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,YON_GT_TSTV_A,4,YON_GT_TSTV_C,
 4,4,4,YON_GT_TSTV_G,4,4,4,4,4,4,4,4,4,4,4,4,YON_GT_TSTV_T,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

struct yon_stats_sample {
public:
	yon_stats_sample(void);
	~yon_stats_sample(void);

	yon_stats_sample& operator+=(const yon_stats_sample& other);

	/**<
	 * Compute the number of transitions and transversions and then use
	 * those numbers to calculate the transition/transversion ratio. These
	 * values are stored internally in the struct.
	 * @return Returns TRUE.
	 */
	bool LazyEvalute(void);

	/**<
	 * Converts this struct into a partial (albeit valid) JSON string.
	 * Takes as argument an input buffer reference and a reference to
	 * the target sample name as described in the global header.
	 * @param buffer      Dst buffer.
	 * @param sample_name Src sample name.
	 * @return            Returns TRUE.
	 */
	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::string& sample_name) const;

	/**<
	 * Resets all internal counters. Useful when reusing
	 * objects without releasing memory.
	 */
	void reset(void);

public:
	uint64_t  n_ins, n_del, n_singleton;
	uint64_t  n_ts, n_tv;
	double    ts_tv_ratio;
	uint64_t* base_conv[9]; // {A,T,G,C,unknown,'.',EOV,ins,del}
	uint64_t* ins_del_dist; // indel distribution. values are zigzagged.
};

struct yon_stats_tstv {
public:
	struct yon_stats_tstv_obj {
		yon_stats_tstv_obj(): n_allele(0), n_non_ref(0), t_non_ref(0), allele_encodings(nullptr), non_ref_encodings(nullptr), b_size(nullptr){}
		yon_stats_tstv_obj(const uint32_t n_allele): n_allele(n_allele), n_non_ref(0), t_non_ref(0), allele_encodings(new uint8_t[n_allele]), non_ref_encodings(new uint8_t[n_allele]), b_size(new int32_t[n_allele]){ }
		~yon_stats_tstv_obj(){
			delete [] allele_encodings;
			delete [] non_ref_encodings;
			delete [] b_size;
		}

		uint32_t  n_allele;
		uint32_t  n_non_ref; // Number of alleles with non-ref alleles.
		uint32_t  t_non_ref; // Target id for non-ref allele.
		uint8_t*  allele_encodings;
		uint8_t*  non_ref_encodings;
		int32_t*  b_size;
	};

public:
	yon_stats_tstv();
	yon_stats_tstv(const uint32_t n_samples);
	~yon_stats_tstv(void);

	yon_stats_tstv& operator+=(const yon_stats_tstv& other);

	// Accessors
	inline yon_stats_sample& operator[](const uint32_t pos){ return(this->sample[pos]); }
	inline const yon_stats_sample& operator[](const uint32_t pos) const{ return(this->sample[pos]); }

	/**<
	 * Addition operator in the situation where the sample order
	 * in the other object is permuted according to a provided
	 * permutation array.
	 * @param other Src object to read data from.
	 * @param ppa   Src permutation array describing the relationship between the current object and the provided one
	 * @return      Returns an reference to the dst object.
	 */
	yon_stats_tstv& Add(const yon_stats_tstv& other, const yon_gt_ppa* ppa);

	/**<
	 * Allocates memory for a number of samples. All previous data
	 * will be deleted without consideration.
	 * @param n_samples
	 */
	void SetSize(const uint32_t n_samples);

	/**<
	 * Construct a valid JSON string from the internal structures. Requires
	 * as input a reference to the global header sample list.
	 * @param buffer       Dst data buffer.
	 * @param sample_names Src vector of sample names.
	 * @return             Returns a reference to the dst data buffer.
	 */
	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::vector<std::string>& sample_names) const;

	/**<
	 * Construct mappings for the alleles into relative offset. Function
	 * will allocate new memory for allele_encodings and non_ref_encodings.
	 * It is not legal to pass the same pointer to allele_encodings and
	 * non_ref_encodings.
	 * @param rcd               Src yon1_t structure that must have genotype data available.
	 * @param allele_encodings  Pointer to empty integer array.
	 * @param non_ref_encodings Pointer to empty integer array.
	 * @return                  Returns TRUE upon success or FALSE otherwise.
	 */
	bool GetEncodings(const yon1_vnt_t& rcd, yon_stats_tstv_obj& helper);

	/**<
	 * Update the current structure with the genotype data from the
	 * provided yon1_t record and a preallocated array of pointers
	 * to yon_gt_rcds. This function requires that the yon1_t record
	 * has genotype data available and that the yon_gt_rcds have been
	 * preprocessed. If genotype data not available or there and/or
	 * there are no meta data available then the function does not
	 * continue. This function is considerably slower than the
	 * Update(const yon1_t&) function that operates directly on the
	 * run-length encoded objects.
	 * @param rcd  Src yon1_t record.
	 * @param rcds Src yon_gt_rcd pointers.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Update(const yon1_vnt_t& rcd, yon_gt_rcd** rcds);

	/**<
	 * Update the current structure with the genotype data from the
	 * provided yon1_t record. This function requires that the yon1_t record
	 * has genotype data available. If genotype data not available or there and/or
	 * there are no meta data available then the function does not
	 * continue. This function does not guarantee that the ordering of the
	 * updates are correct in terms of global ordering. The updates will
	 * be completed according to the (possibly) local permutation array.
	 * Restoring global order is up to the user.
	 * @param rcd  Src yon1_t record.
	 * @param rcds Src yon_gt_rcd pointers.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Update(const yon1_vnt_t& rcd);
	void UpdateDiploid(const yon_gt* gt, yon_stats_tstv_obj& helper);
	void UpdateNPloidy(const yon_gt* gt, yon_stats_tstv_obj& helper);

	// Todo
	void Evaluate(void){
		for(uint32_t i = 0; i < this->n_s; ++i)
			this->sample[i].LazyEvalute();
	}

	/**<
	 * Resets all internal counters. Useful when reusing
	 * objects without releasing memory. Iteratively calls
	 * reset on the child objects.
	 */
	void reset(void);

public:
	uint32_t n_s;
	uint64_t n_rcds;
	uint64_t n_snp, n_mnp, n_ins, n_del, n_other; // determined during lazy evaluation
	uint64_t n_no_alts, n_biallelic, n_multi_allele, n_multi_allele_snp, n_singleton;
	yon_stats_sample* sample;
	uint64_t* alt_count; // number of alt distribution
};

}



#endif /* CORE_TS_TV_OBJECT_H_ */
