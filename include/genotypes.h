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
#ifndef TACHYON_GENOTYPES_H_
#define TACHYON_GENOTYPES_H_

#include <cstring>
#include <unordered_map>

#ifndef SIZE_MAX
#define SIZE_MAX static_cast<size_t>(-1)
#endif

#include "htslib/vcf.h"

#include "buffer.h"
#include "support_vcf.h"
#include "header_footer.h"

namespace tachyon {

/****************************
*  Core genotype
****************************/
#define YON_GT_RCD_MISS 0
#define YON_GT_RCD_EOV  1
#define YON_GT_RCD_REF  2

// Routines for packing/unpacking GT objects.
#define YON_GT_RLE_ALLELE_A(PRIMITITVE, SHIFT, ADD)  (((PRIMITITVE) & ((1 << (SHIFT)) - 1) << (ADD)) >> (ADD));
#define YON_GT_RLE_ALLELE_B(PRIMITIVE, SHIFT, ADD)   (((PRIMITIVE) & ((1 << (SHIFT)) - 1) << ((ADD)+(SHIFT))) >> ((ADD)+(SHIFT)));
#define YON_GT_RLE_LENGTH(PRIMITIVE, SHIFT, ADD)     ((PRIMITIVE) >> (2*(SHIFT) + (ADD)))
#define YON_GT_DIPLOID_ALLELE_LOOKUP(A,B,shift,mask) (((A) & (mask)) << (shift)) | ((B) & (mask))
#define YON_GT_DIPLOID_BCF_A(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> ((SHIFT) + 1)) & (((uint64_t)1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_B(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> 1) & (((uint64_t)1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_PHASE(PRIMITIVE)          ((PRIMITIVE) & 1)
#define YON_GT_BCF1(ALLELE) (((((ALLELE) >> 1) - 1) << 1) | ((ALLELE) & 1))
#define YON_GT_RCD_ALLELE_UNPACK(ALLELE) ((ALLELE) >> 1)
#define YON_GT_RCD_PHASE(ALLELE) ((ALLELE) & 1)

// 0:  Nothing evaluated
// 1:  rcds
// 2:  d_bcf
// 4:  d_bcf_ppa
// 8:  d_exp
// 16: d_occ
#define YON_GT_UN_NONE       0 // nothing
#define YON_GT_UN_RCDS       1 // basic rcds
#define YON_GT_UN_BCF        2 // convert into bcf-style
#define YON_GT_UN_BCF_PPA    4  // convert into bcf-style in-order
#define YON_GT_UN_EXPAND     8 // expand rcds into entries
#define YON_GT_UN_OCC       16 // calculate rcds for each factor
#define YON_GT_UN_ALL       (YON_GT_UN_EXPAND|YON_GT_UN_OCC) // everything

// 0 for missing and 1 for sentinel node. Note that the
// sentinel node never occurs in this encoding type.
const uint8_t YON_GT_RLE_RECODE[3] = {2, 3, 0};

// Vcf:INFO names for fields annotated when triggering
// annotation of genotypes.
const std::vector< std::string > YON_GT_ANNOTATE_FIELDS = {"NM","NPM","AN","HWE_P","AC","AF","AC_P","FS_A","F_PIC","HET","MULTI_ALLELIC"};

/****************************
*  Basic structures
****************************/

/**<
 * Basic structure that maintains the permutation
 * order of the samples in relation to the global header.
 * This object is required if you want to use individual
 * genotypes in the ORIGINAL order. If this is not required
 * in your use-case then this structure has no value.
 */
struct yon_gt_ppa {
public:
	yon_gt_ppa(void);
	yon_gt_ppa(const uint32_t n_samples);
	yon_gt_ppa(const yon_gt_ppa& other);
	yon_gt_ppa(yon_gt_ppa&& other) noexcept;
	yon_gt_ppa& operator=(const yon_gt_ppa& other);
	yon_gt_ppa& operator=(yon_gt_ppa&& other) noexcept;

	~yon_gt_ppa(void);

	uint32_t& operator[](const uint32_t& position){ return(this->ordering[position]); }
	const uint32_t& operator[](const uint32_t& position) const{ return(this->ordering[position]); }
	uint32_t& at(const uint32_t& position){ return(this->ordering[position]); }
	const uint32_t& at(const uint32_t& position) const{ return(this->ordering[position]); }

	/**<
	 * Allocate memory for a new number of samples. All previous
	 * data is deleted without consideration.
	 * @param n_samples Number of samples.
	 */
	void Allocate(const uint32_t n_samples);

	/**<
	 * Restores ordering to [0,..,n_s-1] without clearing data.
	 */
	void reset(void);

	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_gt_ppa& ppa);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_gt_ppa& ppa);

public:
	uint32_t  n_s;
	uint32_t* ordering;
};

/**<
 * Supportive structure in genotype sorter (gtPBWT). Has no
 * other uses.
 */
struct yon_radix_gt {
public:
	yon_radix_gt();
	~yon_radix_gt();
	yon_radix_gt(const yon_radix_gt& other);
	yon_radix_gt(yon_radix_gt&& other);
	yon_radix_gt& operator=(const yon_radix_gt& other);
	yon_radix_gt& operator=(yon_radix_gt&& other);

	bool operator<(const yon_radix_gt& other) const;
	bool operator==(const yon_radix_gt& other) const;
	inline bool operator!=(const yon_radix_gt& other) const{ return(!(*this == other)); }

	friend std::ostream& operator<<(std::ostream& stream, const yon_radix_gt& genotype);

	/**<
	 * Return a bit-packed representation of this data. This should
	 * be considered a private function as it has no external functionality.
	 * @param shift_size Number of bits to shift in.
	 * @return           Returns a bit-packed integer.
	 */
	uint64_t GetPackedInteger(const uint8_t& shift_size = 8) const;

	void resize(const uint8_t new_ploidy);

public:
	uint8_t   n_ploidy;
	uint8_t   n_allocated;
	uint64_t  id;
	uint16_t* alleles;
};

// Primary generic Tachyon FORMAT:GT structure.
/**<
 * Primary Tachyon Format:Gt structure. Encodes alleles
 * as:
 *    0: missing
 *    1: sentinel node (EOV)
 *    2: reference
 *    3: first alt
 *    4: second alt
 *    N: ...
 */
struct yon_gt_rcd {
public:
	yon_gt_rcd();
	yon_gt_rcd(const uint32_t rl, const uint32_t ploidy, const uint8_t* ref) : run_length(rl), allele(new uint8_t[ploidy]){ memcpy(allele, ref, ploidy); }
	yon_gt_rcd(const yon_gt_rcd& other) = delete; // disallow copy ctor (use Clone() instead)
	yon_gt_rcd& operator=(const yon_gt_rcd& other) = delete; // disallow move copy (use Clone() instead)
	yon_gt_rcd(yon_gt_rcd&& other);
	yon_gt_rcd& operator=(yon_gt_rcd&& other);
	~yon_gt_rcd();

	/**<
	 * Cloning subroutine for copying data with the help of externally
	 * provided information. This function requires the knowledge of what
	 * base ploidy this data was constructed with.
	 * @param ploidy Src base ploidy this data was constructed with.
	 * @return       Returns an new instance of yon_gt_rcd.
	 */
	inline yon_gt_rcd Clone(const uint32_t ploidy){ return(yon_gt_rcd(run_length, ploidy, allele)); }

	/**<
	 * Convert this Genotype representation into a valid Vcf string.
	 * @param buffer   Dst buffer.
	 * @param n_ploidy Src base ploidy this data was constructed with.
	 * @return         Returns a reference to the dst buffer.
	 */
	yon_buffer_t& PrintVcf(yon_buffer_t& buffer, const uint8_t& n_ploidy);

public:
	uint32_t run_length;
	uint8_t* allele; // contains phase at first bit
};

// Forward declare.
struct yon_gt_summary;

/**<
 * Genotype container in Tachyon. Stores all genotype records and
 * provides functionality to query and interact with them.
 */
struct yon_gt {
public:
	typedef yonRawIterator<yon_gt_rcd>       iterator;
	typedef yonRawIterator<const yon_gt_rcd> const_iterator;

public:
	yon_gt();
	yon_gt(const yon_gt& other);
	yon_gt& operator=(const yon_gt& other);
	yon_gt(yon_gt&& other) noexcept;
	yon_gt& operator=(yon_gt&& other) noexcept;
	~yon_gt();

	// Iterator
	inline iterator begin(){ return iterator(&this->rcds[0]); }
	inline iterator end()  { return iterator(&this->rcds[this->n_i]); }
	inline const_iterator begin()  const{ return const_iterator(&this->rcds[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->rcds[this->n_i]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->rcds[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->rcds[this->n_i]); }

	/**<
	* Lazy-evaluate internal objects (yon_gt_rcd) from the provided data
	* and method. Converts primitive encoded data into a record structure
	* tuple with (run length, alleles) with the alleles encoded as the first
	* bit for phasing and the remainder seven bits for allele encodings.
	*
	* Decodings are invoked from the EvaluateRecordsM*() functions. These
	* corresponds to the following genotype encodings:
	*     M1) Diploid bi-allelic and no EOV
	*     M2) Diploid any allele count, missingness, and EOV
	*     M4) Nploid any
	*
	* @return Returns TRUE upon success or FALSE otherwise
	*/
	bool Evaluate(void);
	bool EvaluateRecordsM1();
	bool EvaluateRecordsM2();
	bool EvaluateRecordsM4();
	template <class T> bool EvaluateRecordsM1_();
	template <class T> bool EvaluateRecordsM2_();
	template <class T> bool EvaluateRecordsM4_();

	/**<
	 * Wrapper for expanding (possibly) run-length encoded rcds structures
	 * into a vector of length N where each entry corresponds to a single
	 * individual.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool Expand(void);

	/**<
	 * Wrapper for expanding (possibly) run-length encoded rcds structures
	 * into a vector of length N using external memory. No checks are made
	 * to ascertain that the array length is properly allocated. Failure to
	 * allocate sufficient memory will result in a seg-fault.
	 * @param d_expe Dst array of pre-allocated array of length N.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool ExpandExternal(yon_gt_rcd* d_expe);
	bool ExpandRecordsPpa(void);
	bool ExpandRecordsPpaExternal(yon_gt_rcd* d_expe);
	bool ExpandRecords(void);
	bool ExpandRecordsExternal(yon_gt_rcd* d_expe);

	/**<
	* Transform lazy-evaluated Tachyon genotype encodings (d_exp)
	* to htslib bcf1_t genotype encodings.
	* @param rec Input bcf1_t record.
	* @param hdr Input htslib bcf header.
	* @return    Returns the pointer to the input bcf1_t record.
	*/
	bcf1_t* UpdateHtslibGenotypes(bcf1_t* rec, bcf_hdr_t* hdr);

public:
	uint16_t eval_cont;
	uint8_t  add : 7,
             global_phase : 1;
	uint8_t  shift;
	uint8_t  p, m, method; // bytes per entry, base ploidy, base method
	uint32_t n_s, n_i, n_o;     // number samples, number of entries
	uint8_t  n_allele;
	yon_gt_ppa* ppa; // pointer to ppa
	const uint8_t* data; // pointer to data
	yon_gt_rcd* d_exp; // lazy evaluated from ppa/normal to internal offset (length = n_samples). This can be very expensive if evaluated internally for every record.
	yon_gt_rcd* rcds; // lazy interpreted internal records
	uint32_t* n_i_occ;
	yon_gt_rcd** d_occ; // lazy evaluation of occ table
	bool dirty;
};

template <class T>
bool yon_gt::EvaluateRecordsM1_(){
	// Prevent double evaluation.
	if(this->eval_cont & YON_GT_UN_RCDS)
		return true;

	if(this->rcds != nullptr) delete [] this->rcds;
	assert(this->m == 2);
	assert(this->n_allele == 2);

	// Allocate memory for new records.
	this->rcds = new yon_gt_rcd[this->n_i];

	// Reinterpret byte stream into the appropriate actual
	// primitive type.
	const T* r_data = reinterpret_cast<const T*>(this->data);

	// Keep track of the cumulative number of genotypes observed
	// as a means of asserting correctness.
	uint64_t n_total = 0;

	// Iterate over the internal run-length encoded genotypes
	// and populate the rcds structure.
	for(uint32_t i = 0; i < this->n_i; ++i){
		uint8_t phasing = 0;
		if(add) phasing = r_data[i] & 1;
		else    phasing = this->global_phase;

		this->rcds[i].run_length = YON_GT_RLE_LENGTH(r_data[i], shift, add);
		this->rcds[i].allele = new uint8_t[2];
		this->rcds[i].allele[0] = YON_GT_RLE_ALLELE_B(r_data[i], shift, add);
		this->rcds[i].allele[1] = YON_GT_RLE_ALLELE_A(r_data[i], shift, add);
		// Store an allele encoded as (ALLELE << 1 | phasing).
		this->rcds[i].allele[0] = (YON_GT_RLE_RECODE[this->rcds[i].allele[0]] << 1) | phasing;
		this->rcds[i].allele[1] = (YON_GT_RLE_RECODE[this->rcds[i].allele[1]] << 1) | phasing;
		n_total += this->rcds[i].run_length;
	}
	assert(n_total == this->n_s);
	this->eval_cont |= YON_GT_UN_RCDS;
}

template <class T>
bool yon_gt::EvaluateRecordsM2_(){
	// Prevent double evaluation.
	if(this->eval_cont & YON_GT_UN_RCDS)
		return true;

	if(this->rcds != nullptr) delete [] this->rcds;
	assert(this->m == 2);

	// Allocate memory for new records.
	this->rcds = new yon_gt_rcd[this->n_i];

	// Reinterpret byte stream into the appropriate actual
	// primitive type.
	const T* r_data = reinterpret_cast<const T*>(this->data);

	// Keep track of the cumulative number of genotypes observed
	// as a means of asserting correctness.
	uint64_t n_total = 0;

	// Iterate over the internal run-length encoded genotypes
	// and populate the rcds structure.
	for(uint32_t i = 0; i < this->n_i; ++i){
		uint8_t phasing = 0;
		if(add) phasing = r_data[i] & 1;
		else    phasing = this->global_phase;

		this->rcds[i].run_length = YON_GT_RLE_LENGTH(r_data[i], shift, add);
		this->rcds[i].allele = new uint8_t[2];
		this->rcds[i].allele[0] = YON_GT_RLE_ALLELE_B(r_data[i], shift, add);
		this->rcds[i].allele[1] = YON_GT_RLE_ALLELE_A(r_data[i], shift, add);
		// Store an allele encoded as (ALLELE << 1 | phasing).
		this->rcds[i].allele[0] = (this->rcds[i].allele[0] << 1) | phasing;
		this->rcds[i].allele[1] = (this->rcds[i].allele[1] << 1) | phasing;
		n_total += this->rcds[i].run_length;
	}
	assert(n_total == this->n_s);
	this->eval_cont |= YON_GT_UN_RCDS;
}

template <class T>
bool yon_gt::EvaluateRecordsM4_(){
	// Prevent double evaluation.
	if(this->eval_cont & YON_GT_UN_RCDS)
		return true;

	if(this->rcds != nullptr) delete [] this->rcds;
	assert(this->m != 2);

	// Allocate memory for new records.
	this->rcds = new yon_gt_rcd[this->n_i];

	// Keep track of the cumulative number of genotypes observed
	// as a means of asserting correctness.
	uint64_t n_total  = 0;
	uint64_t b_offset = 0;

	// Iterate over the internal run-length encoded genotypes
	// and populate the rcds structure.
	for(uint32_t i = 0; i < this->n_i; ++i){
		const T* run_length = reinterpret_cast<const T*>(&this->data[b_offset]);
		b_offset += sizeof(T);

		this->rcds[i].run_length = *run_length;
		this->rcds[i].allele = new uint8_t[this->m];
		for(uint32_t j = 0; j < this->m; ++j, ++b_offset)
			this->rcds[i].allele[j] = this->data[b_offset];

		n_total += this->rcds[i].run_length;
	}
	assert(n_total == this->n_s);
	this->eval_cont |= YON_GT_UN_RCDS;
}

/****************************
*  Genotype summary statistics
****************************/

/**<
 * Supportive structure used in yon_gt_summary below. Corresponds to the
 * lazy evaluation of a genotype summary statistics ojbect. Pre-computes
 * a variety of standard genotype statistics such as:
 *    1) Base ploidy
 *    2) Allele counts
 *    3) Allele frequencies
 *    4) Number of missing, number of special sentinel symbols
 *    5) Fisher's exact test for strand bias
 *    6) Hardy-Weinberg equilibrium P-value
 *    7) F-statistics for population inbreeding
 *    8) Site heterozygosity
 *
 * This object should be considered internal. Manipulation of this object
 * without modifying all othe related genotype objects could result in
 * significant discord.
 */
struct yon_gt_summary_rcd {
public:
	yon_gt_summary_rcd();
	yon_gt_summary_rcd(const yon_gt_summary_rcd& other);
	yon_gt_summary_rcd& operator=(const yon_gt_summary_rcd& other);
	yon_gt_summary_rcd(yon_gt_summary_rcd&& other) noexcept;
	yon_gt_summary_rcd& operator=(yon_gt_summary_rcd&& other) noexcept;
	~yon_gt_summary_rcd();

public:
	uint8_t n_ploidy;
	uint32_t n_ac_af;
	uint32_t n_fs;
	uint32_t* ac; // allele counts
	float* af; // allele frequency
	uint32_t nm; // number of non-sentinel, non-missing symbols
	uint32_t npm;// number of missing symbols
	uint32_t an; // number of non-sentinel symbols
	float* fs_a;// fisher strand test p
	float hwe_p;// hardy-weinberg p
	float f_pic;
	float heterozygosity;
};

/**<
 * Primary object used in the genotype trie in yon_gt_summary below. This object is built
 * using information about the number of alleles. However this information
 * is not explicitly stored here (it is stored in the parent container).
 * Because of this it is not possible to invoke the copy ctor or copy assign
 * operator (use Clone() instead). It also prevents the use of this object
 * in direct arithmetic use cases.
 */
struct yon_gt_summary_obj {
public:
	yon_gt_summary_obj();
	yon_gt_summary_obj(const yon_gt_summary_obj& other) = delete;
	yon_gt_summary_obj& operator=(const yon_gt_summary_obj& other) = delete;
	yon_gt_summary_obj(yon_gt_summary_obj&& other) noexcept;
	yon_gt_summary_obj& operator=(yon_gt_summary_obj&& other) noexcept;
	~yon_gt_summary_obj();

	// Clone helper ctor.
	yon_gt_summary_obj(const uint8_t n_alleles, const uint64_t cnt, yon_gt_summary_obj* c);

	yon_gt_summary_obj Clone(const uint8_t n_alleles){ return(yon_gt_summary_obj(n_alleles, n_cnt, children)); }

	inline yon_gt_summary_obj& operator[](const uint32_t pos){ return(this->children[pos]); }
	inline const yon_gt_summary_obj& operator[](const uint32_t pos) const{ return(this->children[pos]); }

public:
	uint64_t n_cnt;
	yon_gt_summary_obj* children;
};

struct yon_gt_summary {
public:
	yon_gt_summary(void);
	yon_gt_summary(const uint8_t base_ploidy, const uint8_t n_alleles);
	yon_gt_summary(yon_gt_summary&& other) noexcept;
	yon_gt_summary& operator=(yon_gt_summary&& other) noexcept;
	yon_gt_summary(const yon_gt_summary& other);
	yon_gt_summary& operator=(const yon_gt_summary& other);
	~yon_gt_summary();

	void Setup(const uint8_t base_ploidy, const uint8_t n_als);

	/**<
	 * Iteratively add layers to the complete trie in order to represent
	 * a possible ploidy-dimensional matrix at the leafs. The memory
	 * cost of the trie is O(n_alleles ^ ploidy) and allows the prefix
	 * lookup of any partial genotype.
	 *children
	 * @param target
	 * @param depth
	 * @return
	 */
	bool AddGenotypeLayer(yon_gt_summary_obj* target, uint8_t depth);

	yon_gt_summary& operator+=(const yon_gt& gt);
	yon_gt_summary& Add(const yon_gt& gt, const uint32_t n_i, const yon_gt_rcd* rcds);

	// Accessors to internal data.
	inline uint32_t* GetAlleleCountsRaw(void){ return(this->alleles); }
	inline const uint32_t* GetAlleleCountsRaw(void) const{ return(this->alleles); }

	/**<
	 * Calculates and returns a vector of alelle counts.
	 * @return Returns a vector of allele counts.
	 */
	std::vector<uint64_t> GetAlleleCounts(void) const;

	/**<
	 * Calculates and returns a vector of tuples (allele count, allele frequency)
	 * for each unique allele.
	 * @return Returns a vector of tuples (allele count, allele frequency) of alleles.
	 */
	std::vector< std::pair<uint64_t,double> > GetAlleleCountFrequency(void) const;

	/**<
	 * Calculates and returns a vector of vectors of allele counts. The structure is
	 * x[ploidy][allele] corresponding to, for example, maternal (0) and paternal (1)
	 * chromosome for a given allele.
	 * @return Returns a vector of vectors of allele counts of allle counts.
	 */
	std::vector< std::vector<uint64_t> > GetAlleleStrandCounts(void) const;

	/**<
	 * Traverse the genotypic trie and return the
	 * @param data
	 * @param target
	 * @param depth
	 * @return
	 */
	bool GetGenotype(std::vector<uint64_t>& data,
	                 yon_gt_summary_obj* target,
	                 uint8_t depth) const;

	// Todo: unfinished.
	std::vector<yon_gt_rcd> GetGenotypeCounts(bool drop_empty = true) const;

	/**<
	 * Calculates and returns a vector of Fisher's exact test P-values for strand
	 * bias. This function requires that the data is diploid.
	 * @param phred_scale Predicate for PHRED-scaling the data.
	 * @return            Returns a vector of Fisher's exact P-values or empty if not diploid.
	 */
	std::vector<double> GetStrandBiasAlleles(const bool phred_scale = true) const;

	/**<
	 * Calulcates and returns the Hardy-Weinberg P-value. This function requires the
	 * data to be diploid and biallelic.
	 * @return Returns the Hardy-Weinberg P-value or -1 if not diploid and biallelic.
	 */
	double CalculateHardyWeinberg(void) const;

	/**<
	 * Lazy evaluates the current data in this object into easy-to-use
	 * elements of the yon_gt_summary_rcd structure.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool LazyEvaluate(void);

public:
	uint8_t    n_ploidy; // base ploidy at site
	uint8_t    n_alleles; // number of alleles
	uint32_t*  alleles; // allele counts
	uint32_t** alleles_strand; // allelic counts per chromosome
	yon_gt_summary_obj* gt; // complete genotypic trie with branch-size n_alleles
	yon_gt_summary_rcd* d; // lazy evaluated record
};

/****************************
*  Occ functionality
****************************/

/**<
 * Structure for using the partial-sum algortihms with constrained run-length
 * encoded genotypes. Allows for O(1)-time partitions into one or more groupings.
 */
struct yon_occ {
public:
	typedef std::unordered_map<std::string, uint32_t> map_type;

public:
	yon_occ() = default;
	~yon_occ() = default;

	/**<
	 * Read a target file from disk and parse its contents into a valid Occ
	 * table.
	 * @param file_name Src file-name to read.
	 * @param header    Src variant header reference.
	 * @param delimiter Delimiter character used in the src file.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool ReadTable(const std::string file_name, const yon_vnt_hdr_t& header, const char delimiter = '\t');

	/**<
	 * Construct an Occ table from the pre-loaded matrix of sample->groupings.
	 * The BuildTable() and BuildTable(yon_gt_ppa*) function requires that
	 * ReadTable() has been run and successfully completed.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool BuildTable(void);
	bool BuildTable(const yon_gt_ppa* ppa_p);

public:
	// Map from group name to row offset in the table.
	map_type map;
	// Unique names of grouping factors.
	std::vector<std::string> row_names;
	// Total cumulative sums for each row.
	std::vector<uint32_t> cum_sums;

	// A matrix with proportions samples times groupings
	// rows corresponds to the cumulative sum of a grouping
	// over the samples. The table corresponds to the set
	// membership (presence or absence) and the occ table
	// corresponds to the cumsum at any given sample offset.
	std::vector< std::vector<uint32_t> > table;
	std::vector< std::vector<uint32_t> > occ;
	// Transpose of occ table for data locality lookups when
	// using constrained run-length encoded objects.
	std::vector< std::vector<uint32_t> > vocc;
};

/****************************
*  Supportive structures
****************************/

/**<
 * Primary bit-packed controller for a yon1_vnt_t record. Keeps track
 * of a variety of essential boolean fields required for the proper
 * functioning of tachyon.
 */
struct yon_vnt_cnt {
public:
	yon_vnt_cnt(void);
	~yon_vnt_cnt() = default;

	friend yon_buffer_t& operator+=(yon_buffer_t& buffer, const yon_vnt_cnt& entry){
		buffer += (uint16_t)*reinterpret_cast<const uint16_t* const>(&entry);
		return(buffer);
	}

	/**<
	 * Interpret a packed uint16_t into a controller struct and
	 * overload this object with it.
	 * @param value Src packed integer.
	 */
	void operator=(const uint16_t& value);

	/**<
	 * Convert this structure into a bit-packed uint16_t value.
	 * @return Returns a bit-packed representation of this object.
	 */
	inline uint16_t ToValue(void) const{ return((uint16_t)*reinterpret_cast<const uint16_t* const>(this)); }

public:
	/**< Controller field. The first seven fields describes
	 * genotype-specific information. The remainder bit-fields
	 * describes variant-specific information.
	 * 1) If genotype data is available
	 * 2) If there is are any missing genotypes
	 * 3) The phase of all genotypes if diploid and no mixed phase
	 * 4) If all genotypes share the same phase
	 * 5) What compression method genotypes are stored in
	 * 6) What primitive type the genotypes are encoded as
	 * 7) If this variant is biallelic
	 * 8) If this variant is a simple SNP->SNP
	 * 9) If the genotypes are base diploid
	 * 10) If there is mixed ploidy in the genotypes (i.e. EOV values present)
	 * 11) If allele data is bit-packed
	 * 12) If all alleles are SNVs
	 */
	uint16_t
		gt_available:        1, // if there is any GT data
		gt_has_missing:      1, // any missing
		gt_phase_uniform:    1, // all phased/unphased
		gt_has_mixed_phasing:1, // has mixed phasing
		gt_compression_type: 4, // GT compression algorithm
		gt_primtive_type:    2, // type of compression primitive (unsinged integer types)
		gt_mixed_ploidy:     1, // has mixed ploidy (have a sentinel node symbol)
		biallelic:           1, // is biallelic
		simple_snv:          1, // is simple SNV->SNV
		diploid:             1, // is diploid
		alleles_packed:      1, // are the alleles packed in a BYTE
		all_snv:             1; // are all ref/alt simple SNVs
};

/****************************
*  Supportive conversion functions
****************************/
template <class T>
static yon_gt* GetGenotypeDiploidRLE(
	const char* data,
	const uint32_t n_entries,
	const uint32_t n_samples,
	const uint16_t n_als,
	const uint8_t  n_ploidy,
	const uint16_t controller,
	yon_gt_ppa* ppa)
{
	yon_gt* x = new yon_gt;
	const yon_vnt_cnt* cont = reinterpret_cast<const yon_vnt_cnt*>(&controller);
	x->shift  = cont->gt_has_missing     ? 2 : 1;
	x->add    = cont->gt_has_mixed_phasing ? 1 : 0;
	x->global_phase = cont->gt_phase_uniform;
	x->data   = reinterpret_cast<const uint8_t*>(data);
	x->n_i    = n_entries;
	x->method = 1;
	x->p     = sizeof(T);
	x->m     = 2;
	x->n_s   = n_samples;
	x->n_allele = n_als;
	x->ppa = ppa;
	return(x);
}

template <class T>
static yon_gt* GetGenotypeDiploidSimple(
	const char* data,
	const uint32_t n_entries,
	const uint32_t n_samples,
	const uint16_t n_als,
	const uint8_t  n_ploidy,
	const uint16_t controller,
	yon_gt_ppa* ppa)
{
	yon_gt* x = new yon_gt;
	const yon_vnt_cnt* cont = reinterpret_cast<const yon_vnt_cnt*>(&controller);
	x->n_allele = n_als;
	x->shift = ceil(log2(x->n_allele + 2 + 1));
	x->add   = cont->gt_has_mixed_phasing ? 1 : 0;
	x->global_phase = cont->gt_phase_uniform;
	x->data = reinterpret_cast<const uint8_t*>(data);
	x->n_i = n_entries;
	x->method = 2;
	x->p = sizeof(T);
	x->m = 2;
	x->n_s = n_samples;
	x->ppa = ppa;

	return(x);
}

template <class T>
static yon_gt* GetGenotypeNploid(
	const char* data,
	const uint32_t n_entries,
	const uint32_t n_samples,
	const uint16_t n_als,
	const uint8_t  n_ploidy,
	const uint16_t controller,
	yon_gt_ppa* ppa)
{
	yon_gt* x = new yon_gt;
	const yon_vnt_cnt* cont = reinterpret_cast<const yon_vnt_cnt*>(&controller);
	x->shift = 0;
	x->add   = cont->gt_has_mixed_phasing;
	x->global_phase = cont->gt_phase_uniform;
	x->data = reinterpret_cast<const uint8_t*>(data);
	x->n_i = n_entries;
	x->method = 4;
	x->m = n_ploidy;
	x->p = sizeof(T);
	x->n_s = n_samples;
	x->n_allele = n_als;
	x->ppa = ppa;

	return(x);
}

}

#endif
