/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_VARIANT_RECORD_H_
#define TACHYON_VARIANT_RECORD_H_

#include "primitive_container.h"
#include "genotypes.h"

#include <cmath>

namespace tachyon {

/****************************
*  Core
****************************/
const char* const YON_REFALT_LOOKUP = "ATGC.XN";
#define YON_ALLELE_A       0
#define YON_ALLELE_T       1
#define YON_ALLELE_G       2
#define YON_ALLELE_C       3
#define YON_ALLELE_MISS    4
#define YON_ALLELE_NON_REF 5
#define YON_ALLELE_N       6

/****************************
*  TS/TV objects
****************************/
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


struct yon_allele {
public:
	yon_allele();
	yon_allele(const std::string& string);
	yon_allele& operator=(const std::string& value);
	yon_allele& operator=(const char* value);
	yon_allele& operator=(const char& value);
	yon_allele(const yon_allele& other);
	yon_allele(yon_allele&& other) noexcept;
	yon_allele& operator=(const yon_allele& other);
	yon_allele& operator=(yon_allele&& other) noexcept;
	~yon_allele();

	yon_allele& ParseFromBuffer(const char* const in);

	inline const uint16_t& size(void) const{ return(this->l_allele); }
	inline const uint16_t& length(void) const{ return(this->l_allele); }
	inline const std::string ToString(void) const{ return(std::string(this->allele, this->l_allele)); }

	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_allele& entry);

public:
	uint16_t l_allele;
	char*    allele;
};

/**<
 * Primary evaluated record of a Tachyon variant. Construction is done
 * outside of this definition. Evaluation into this object is relatively
 * inexpensive in isolation but quiet expensive when considered
 * across the an entire dataset as this structure needs be evaluated many
 * millions of times over many different containers.
 */
struct yon1_vnt_t {
public:
	yon1_vnt_t();
	yon1_vnt_t(const yon1_vnt_t& other);
	yon1_vnt_t& operator=(const yon1_vnt_t& other);
	yon1_vnt_t(yon1_vnt_t&& other) noexcept;
	//yon1_vnt_t& operator=(yon1_vnt_t&& other) noexcept = delete; // temporarily deleted
	~yon1_vnt_t();

	bool UpdateBase(const bcf1_t* record);

	bool EvaluateSummary(bool lazy_evaluate = true);
	bool EvaluateOcc(yon_occ& occ);
	bool EvaluateOccSummary(bool lazy_evaluate = true);

	/**<
	 * Check if it is possible to pack the REF and ALT allele strings into
	 * a single uint8_t. The input allelic data has to be diploid, biallelic and
	 * match the regular expression pattern "^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$"
	 * @return Returns TRUE if it is possible to pack the ref and alt alleles into a byte or FALSE otherwise.
	 */
	bool UsePackedRefAlt(void) const;

	/**<
	 * Bitpack biallelic, diploid REF and ALT data into a single uint8_t. Failure
	 * to check for validity beforehand with UsePackedRefAlt() may result in
	 * errors.
	 * @return Returns a bitpacked unsigned byte.
	 */
	uint8_t PackRefAltByte(void) const;

	/**<
	 * Returns the alleles as a concatenated string. For example the alleles
	 * A and T is returned as "A,T". This function is used primarily in the
	 * UpdateHtslibVcfRecord() function for exporting into a htslib bcf1_t record.
	 * @return Returns a concatenated string of alleles.
	 */
	std::string GetAlleleString(void) const;

	// Base setters for function pointer overloading.
	inline void SetController(const uint16_t value){ this->controller = value; }
	inline void SetBasePloidy(const uint8_t value){ this->n_base_ploidy = value; }
	inline void SetInfoPatternId(const int32_t value){ this->info_pid = value; }
	inline void SetFormatPatternId(const int32_t value){ this->fmt_pid = value; }
	inline void SetFilterPatternId(const int32_t value){ this->flt_pid = value; }
	inline void SetQuality(const float value){ this->qual = value; }
	inline void SetPosition(const int64_t value){ this->pos = value; }
	inline void SetChromosome(const uint32_t value){ this->rid = value; }
	inline void SetName(const std::string& value);

	/**<
	 * Updates a htslib bcf1_t record with data available in this meta record.
	 * This function is used when converting yon1_t records to bcf1_t records.
	 * @param rec Input bcf1_t record that has been allocated.
	 * @param hdr Input bcf hdr structure converted from tachyon header.
	 * @return Returns the input bcf1_t record pointer.
	 */
	bcf1_t* UpdateHtslibVcfRecord(bcf1_t* rec, bcf_hdr_t* hdr) const;

	void OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr);

	void OutputHtslibVcfFormat(bcf1_t* rec,
	                           bcf_hdr_t* hdr,
	                           DataBlockSettings& settings,
	                           yon_gt_rcd* external_exp) const;

	void OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr) const;

	void ToVcfString(const VariantHeader& header,
	           yon_buffer_t& buffer,
	           uint32_t display,
	           yon_gt_rcd* external_rcd = nullptr) const;

	bool AddInfoFlag(const std::string& tag, VariantHeader& header);

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, int_t& data_point){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = new PrimitiveContainer<int_t>(data_point);
		} else {
			info[n_info++] = new PrimitiveContainer<int_t>(data_point);
			info_hdr.push_back(info_tag);
		}
		return true;
	}

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, int_t* data_point, const size_t n_points){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = new PrimitiveContainer<int_t>(data_point, n_points);
		} else {
			info[n_info++] = new PrimitiveContainer<int_t>(data_point, n_points);
			info_hdr.push_back(info_tag);
		}
		return true;
	}

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, PrimitiveContainer<int_t>*& container){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = container;
		} else {
			info[n_info++] = container;
			info_hdr.push_back(info_tag);
		}
		container = nullptr;
		return true;
	}

	bool AddFilter(const std::string& tag, VariantHeader& header);
	bool AddGenotypeStatistics(VariantHeader& header, const bool replace_existing = true);
	bool AddGenotypeStatisticsOcc(VariantHeader& header, std::vector<std::string>& names);

	/**<
	 * Calculates various F-statistics over the provided groupings using
	 * the Occ table. Calculates one set of statistics over all the groupings
	 * and for each individual group. The calculations here are available only
	 * for diploid individuals at biallelic sites.
	 *
	 * Warning: this function returns incorrect results if the provided
	 *          grouping sets are overlapping. The incorrect results arise
	 *          from the fact that a single individual cannot simultaneously
	 *          exist in multiple distinct populations.
	 * @param header Src VariantHeader.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool EvaluatePopGenOcc(VariantHeader& header);

	PrimitiveContainerInterface* GetInfo(const std::string& name) const;
	PrimitiveGroupContainerInterface* GetFmt(const std::string& name) const;

	const VcfFilter* GetFlt(const std::string& name) const;
	int32_t GetInfoOffset(const std::string& name) const;
	int32_t GetFormatOffset(const std::string& name) const;
	int32_t GetFilterOffset(const std::string& name) const;

public:
	bool is_dirty; // if data has been modified
	bool is_loaded_gt;
	yon_vnt_cnt controller;
	uint8_t  n_base_ploidy;
	uint16_t n_alleles;
	int32_t  n_flt, n_fmt, n_info;
	int32_t  info_pid, flt_pid, fmt_pid;
	int32_t  m_fmt, m_info, m_allele; // allocated sizes
	float    qual;
	uint32_t rid;
	int64_t  pos;
	std::string name;

	std::vector<const YonInfo*> info_hdr;
	std::vector<const YonFormat*> fmt_hdr;
	std::vector<const VcfFilter*> flt_hdr;

	std::unordered_map<std::string, uint32_t> info_map;
	std::unordered_map<std::string, uint32_t> fmt_map;
	std::unordered_map<std::string, uint32_t> flt_map;

	yon_allele* alleles;
	yon_gt* gt;
	yon_gt_summary* gt_sum;
	yon_gt_summary* gt_sum_occ; // summary if occ has been invoked
	PrimitiveContainerInterface** info;
	PrimitiveGroupContainerInterface** fmt;
};

// Genotype helper
struct GenotypeSummary {
public:
	GenotypeSummary(void) :
		base_ploidy(0),
		phase_if_uniform(0),
		mixed_phasing(false),
		invariant(false),
		n_missing(0),
		n_vector_end(0)
	{}

	~GenotypeSummary() = default;

	/**<
	 * Gathers summary statistics for a vector of genotypes
	 * at a given site. Collects information regarding the
	 * number of missing genotypes and count of sentinel
	 * nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 * Todo: Use hashes to check for uniformity of genotypes.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @param fmt       The target htslib format structure.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	template<class T>
	bool Evaluate(const size_t& n_samples, const bcf_fmt_t& fmt){
		if(fmt.p_len == 0) return true;
		assert(fmt.size/fmt.n == sizeof(T));

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = fmt.n;

		// Find first phase
		this->phase_if_uniform = 0;
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		int j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(io::VcfGenotype<T>::IsMissing(fmt.p[j]) == true
			   || io::VcfType<T>::IsVectorEnd(fmt.p[j]) == true)
				j += fmt.n;
			else {
				this->phase_if_uniform = fmt.p[j] & 1;
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		// symbols and assess uniformity of phasing.
		j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(io::VcfGenotype<T>::IsMissing(fmt.p[j]) == false
			   && io::VcfType<int8_t>::IsVectorEnd(fmt.p[j]) == false
			   && (fmt.p[j] & 1) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < fmt.n; ++k, ++j){
				this->n_missing    += io::VcfGenotype<T>::IsMissing(fmt.p[j]);
				this->n_vector_end += io::VcfType<T>::IsVectorEnd(fmt.p[j]);
			}
		}

		return true;
	}

	/**<
	 * Gathers summary statistics for a vector of genotypes at a given site.
	 * Collects information regarding the number of missing genotypes and
	 * count of sentinel nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 *
	 * Todo: Use hashes to check for uniformity of genotypes.
	 *
	 * @param rec       The src yon1_vnt_t record.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool Evaluate(const yon1_vnt_t& rec, const size_t& n_samples){
		if(rec.gt == nullptr)
			return false;

		if(rec.gt->rcds == nullptr)
			return false;

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = rec.gt->m;

		// Find first phase
		this->phase_if_uniform = 0;
		const uint8_t p_offset = std::max(0, (int)rec.gt->m - 1);
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		for(uint32_t i = 0; i < rec.gt->n_i; ++i){
			if(YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[p_offset]) > 1){
				this->phase_if_uniform = (YON_GT_RCD_PHASE(rec.gt->rcds[i].allele[p_offset]));
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		for(uint32_t i = 0; i < rec.gt->n_i; ++i){
			if(YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[p_offset]) > 1
			   && (YON_GT_RCD_PHASE(rec.gt->rcds[i].allele[p_offset])) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < rec.gt->m; ++k){
				this->n_missing    += YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[k]) == YON_GT_RCD_MISS ? rec.gt->rcds[i].run_length : 0;
				this->n_vector_end += YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[k]) == YON_GT_RCD_EOV ? rec.gt->rcds[i].run_length : 0;
			}
		}

		//std::cerr << utility::timestamp("DEBUG") << n_missing << "," << n_vector_end << "," << (int)base_ploidy << "," << phase_if_uniform << ":" << mixed_phasing << ":" << invariant << std::endl;

		return true;
	}

	inline bool isBaseDiploid(void) const{ return(this->base_ploidy == 2); }

public:
	uint8_t  base_ploidy;
	bool     phase_if_uniform;
	bool     mixed_phasing;
	bool     invariant;
	uint64_t n_missing;
	uint64_t n_vector_end;
};

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
	yon_buffer_t& ToJsonString(yon_buffer_t& buffer, const std::string& sample_name) const;

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

/**<
 * Supportive structure for computing variant-centric summary statistics
 * using the genotypes and site-specific meta information.
 */
struct yon_stats_tstv {
public:
	// Supportive structure.
	struct yon_stats_tstv_obj {
		yon_stats_tstv_obj(): n_allele(0), n_non_ref(0), t_non_ref(0), allele_encodings(nullptr), non_ref_encodings(nullptr), b_size(nullptr){}
		yon_stats_tstv_obj(const uint32_t n_allele): n_allele(n_allele), n_non_ref(0), t_non_ref(0), allele_encodings(new uint8_t[n_allele]), non_ref_encodings(new uint8_t[n_allele]), b_size(new int32_t[n_allele]){ }
		~yon_stats_tstv_obj(){
			delete [] allele_encodings;
			delete [] non_ref_encodings;
			delete [] b_size;
		}

		uint32_t n_allele;
		uint32_t n_non_ref; // Number of alleles with non-ref alleles.
		uint32_t t_non_ref; // Target id for non-ref allele.
		uint8_t* allele_encodings;
		uint8_t* non_ref_encodings;
		int32_t* b_size;
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
	yon_buffer_t& ToJsonString(yon_buffer_t& buffer, const std::vector<std::string>& sample_names) const;

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

	/**<
	 * Short-hand support function for invoking LazyEvaluate() in all available
	 * children objects. See LazyEvaluate() in the children structure for more
	 * information.
	 */
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

#endif /* CORE_VARIANT_RECORD_H_ */
