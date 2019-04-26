#ifndef INDEX_VARIANT_INDEX_H_
#define INDEX_VARIANT_INDEX_H_

#include <cstring>
#include <cmath>
#include <vector>

#include "variant_index_bin.h"
#include "variant_index_contig.h"
#include "variant_index_linear.h"
#include "header_footer.h"
#include "io/vcf_utils.h"
#include "generic_iterator.h"

namespace tachyon {
namespace index {

class VariantIndexQuadTree {
public:
	typedef VariantIndexQuadTree self_type;
    typedef std::size_t          size_type;
    typedef VariantIndexContig   value_type;
    typedef VariantIndexLinear   linear_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef yon1_idx_rec         linear_entry_type;
    typedef YonContig            contig_type;
    typedef linear_type*         linear_pointer;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
   	VariantIndexQuadTree();
   	VariantIndexQuadTree(const self_type& other);
	~VariantIndexQuadTree();

	self_type& operator+=(const self_type& other) {
		assert(this->n_contigs_ == other.n_contigs_);
		for (int i = 0; i < this->n_contigs_; ++i)
			this->contigs_[i] += other.contigs_[i];

		return(*this);
	}

	// Element access
	inline reference at(const size_type& position) { return(this->contigs_[position]); }
	inline const_reference at(const size_type& position) const { return(this->contigs_[position]); }
	inline reference operator[](const size_type& position) { return(this->contigs_[position]); }
	inline const_reference operator[](const size_type& position) const { return(this->contigs_[position]); }
	inline pointer data(void) { return(this->contigs_); }
	inline const_pointer data(void) const { return(this->contigs_); }
	inline reference front(void) { return(this->contigs_[0]); }
	inline const_reference front(void) const { return(this->contigs_[0]); }
	inline reference back(void) { return(this->contigs_[this->n_contigs_ - 1]); }
	inline const_reference back(void) const { return(this->contigs_[this->n_contigs_ - 1]); }

	// Capacity
	inline bool empty(void) const { return(this->n_contigs_ == 0); }
	inline const size_type& size(void) const { return(this->n_contigs_); }
	inline const size_type& capacity(void) const { return(this->n_capacity_); }

	// Iterator
	inline iterator begin() { return iterator(&this->contigs_[0]); }
	inline iterator end() { return iterator(&this->contigs_[this->n_contigs_]); }
	inline const_iterator begin() const { return const_iterator(&this->contigs_[0]); }
	inline const_iterator end() const { return const_iterator(&this->contigs_[this->n_contigs_]); }
	inline const_iterator cbegin() const { return const_iterator(&this->contigs_[0]); }
	inline const_iterator cend() const { return const_iterator(&this->contigs_[this->n_contigs_]); }

	self_type& Add(const std::vector<VcfContig>& contigs);
	self_type& Add(const std::vector<YonContig>& contigs);

	/**<
	 * Add a contig with n_levels to the chain
	 * @param l_contig Length of contig
	 * @param n_levels Number of desired 4^N levels
	 * @return         Returns a reference of self
	 */
	inline self_type& Add(const uint32_t& contigID, const uint64_t& l_contig, const uint8_t& n_levels) {
		if (this->size() + 1 >= this->n_capacity_)
			this->resize();

		new( &this->contigs_[this->n_contigs_] ) value_type( contigID, l_contig, n_levels );
		++this->n_contigs_;
		return(*this);
	}

	/**<
	 * Resizes the index to accept more contigs than currently allocated
	 * memory for. Resizes for the quad-tree index and the linear index
	 */
	void resize(void);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& index);
	friend std::istream& operator>>(std::istream& stream, self_type& index);

public:
	size_type n_contigs_; // number of contigs
	size_type n_capacity_;
	pointer   contigs_;
};

}
}

#endif /* INDEX_VARIANT_INDEX_H_ */
