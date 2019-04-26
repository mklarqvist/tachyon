#ifndef CONTAINERS_VCF_CONTAINER_H_
#define CONTAINERS_VCF_CONTAINER_H_

#include "generic_iterator.h"
#include "io/vcf_utils.h"
#include "variant_record.h"
#include "io/vcf_reader.h"

namespace tachyon {
namespace containers {

/**<
 * This class is both a container and iterator for htslib vcf entries. Primarily
 * used during import from Vcf/Bcf into a Tachyon archive. Retrieves Vcf records
 * from a byte stream or disk if an active VcfReader is provided.
 */
class VcfContainer{
public:
	typedef VcfContainer       self_type;
	typedef bcf1_t             value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
	VcfContainer(void);
	VcfContainer(const size_type& start_capacity);
	VcfContainer(const VcfContainer& other) = delete; // Disallow copy ctor
	VcfContainer(self_type&& other) noexcept;
	~VcfContainer();
	VcfContainer& operator=(self_type&& other) noexcept;
	VcfContainer& operator=(const self_type& other) = delete; // Disallow assign copy

	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline size_type sizeWithoutCarryOver(void) const{ return(this->n_entries_ - this->n_carry_over_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }
	inline pointer front(void) { return(this->entries_[0]); }
	inline const_pointer front(void) const{ return(this->entries_[0]); }
	inline pointer back(void) { return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }
	inline const_pointer back(void) const{ return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }

	inline void operator+=(const pointer entry) { this->entries_[this->n_entries_++] = entry; }
	inline pointer operator[](const uint32_t position) { return(this->entries_[position]); }
	inline const_pointer operator[](const uint32_t position) const{ return(this->entries_[position]); }
	inline pointer at(const uint32_t position) { return(this->entries_[position]); }
	inline const_pointer at(const uint32_t position) const{ return(this->entries_[position]); }

	inline pointer end(void) { return(this->entries_[this->n_entries_]); }
	inline const_pointer end(void) const{ return(this->entries_[this->n_entries_]); }

	void resize(const size_t new_size);

	/**<
	 * Read and, optionally, parse htslib bcf1_t entries using the provided
	 * VcfReader. Continues to load records until reaching either a fixed
	 * number of variants or number of base-pairs, whichever comes first.
	 * Automatically breaks if a record maps to a different contig relative
	 * the first one retrieved. This is required for importing into a sorted
	 * Tachyon archive.
	 * @param n_variants   Maximum number of records to load.
	 * @param n_bases      Maximum number of base-pairs to load.
	 * @param reader       VcfReader object with an active input stream.
	 * @param unpack_level htslib unpack level (see BCF_UN_ALL).
	 * @return             Returns TRUE if any data was loaded or FALSE otherwise.
	 */
	bool GetVariants(const int32_t n_variants,
	                 const int64_t n_bases,
	                 std::unique_ptr<io::VcfReader>& reader,
	                 const int unpack_level = BCF_UN_ALL);

	/**<
	 * Calculate genotype summary statistics from a lazy evaluated bcf1_t struct.
	 * Warning: this function does NOT check if the Format:GT field exists either
	 * in the header or in the structure itself. The assumption is that it does
	 * exist and according to the Bcf specification has to be the first Format
	 * field set.
	 * @param position  Array offset to src htslib bcf1_t entry in this container.
	 * @param n_samples Number of samples as described in the Vcf header.
	 * @return          Returns a VcfGenotypeSummary object.
	 */
	GenotypeSummary GetGenotypeSummary(const uint32_t position, const uint64_t& n_samples) const;

	/**<
	 * Clear data without releasing allocated memory.
	 */
	void clear(void);

public:
	uint32_t  n_carry_over_;
	size_type n_entries_;
	size_type n_capacity_;
	pointer*  entries_;
};

}
}



#endif /* CONTAINERS_VCF_CONTAINER_H_ */
