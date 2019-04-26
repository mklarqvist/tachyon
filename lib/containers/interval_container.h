#ifndef CONTAINERS_INTERVAL_CONTAINER_H_
#define CONTAINERS_INTERVAL_CONTAINER_H_

#include <regex>

#include "tachyon.h"
#include "utility.h"
#include "third_party/intervalTree.h"
#include "index.h"
#include "header_footer.h"
#include "variant_record.h"
#include "generic_iterator.h"

namespace tachyon{
namespace containers{

class IntervalContainer {
public:
	typedef IntervalContainer      self_type;
    typedef std::size_t            size_type;
    typedef algorithm::Interval<uint32_t, int64_t> interval_type;
    typedef algorithm::IntervalTree<uint32_t, int64_t>  value_type;
    typedef value_type&            reference;
    typedef const value_type&      const_reference;
    typedef value_type*            pointer;
    typedef const value_type*      const_pointer;
    typedef yon_index_t            index_type;
    typedef yon1_idx_rec           index_entry_type;
    typedef yon_vnt_hdr_t          header_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    IntervalContainer();
    IntervalContainer(const self_type& other);
    IntervalContainer(self_type&& other) noexcept;
    ~IntervalContainer(void);
    IntervalContainer& operator=(const self_type& other);
    IntervalContainer& operator=(self_type&& other) noexcept;

    // Element access
    inline reference at(const size_type& position) { return(this->__entries[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
    inline reference operator[](const size_type& position) { return(this->__entries[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
    inline pointer data(void) { return(this->__entries); }
    inline const_pointer data(void) const{ return(this->__entries); }
    inline reference front(void) { return(this->__entries[0]); }
    inline const_reference front(void) const{ return(this->__entries[0]); }
    inline reference back(void) { return(this->__entries[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries_ - 1]); }

    // Capacity
    inline bool empty(void) const{ return(this->n_entries_ == 0); }
    inline const size_type& size(void) const{ return(this->n_entries_); }

    // Iterator
    inline iterator begin() { return iterator(&this->__entries[0]); }
    inline iterator end() { return iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries_]); }

	inline bool HasBlockList(void) const{ return(this->block_list_.size()); }
	inline std::vector<index_entry_type>& GetBlockList(void) { return(this->block_list_); }
	inline const std::vector<index_entry_type>& GetBlockList(void) const{ return(this->block_list_); }
	inline bool HasIntervals(void) const{ return(this->interval_list_.size()); }

    // Interpret
    bool ValidateIntervalStrings(std::vector<std::string>& interval_strings);

    /**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool ParseIntervals(std::vector<std::string>& interval_strings,
	                    const header_type& header,
	                    const index_type& index);

	bool Build(const header_type& header);

	inline std::vector<interval_type> FindOverlaps(const uint32_t& contigID,
	                                               const int64_t& start_position,
	                                               const int64_t& end_position) const
	{
		if (contigID > this->size()) return(std::vector<interval_type>());
		return(this->at(contigID).findOverlapping(start_position, end_position));
	}

	inline std::vector<interval_type> FindOverlaps(const yon1_vnt_t& rcd) const {
		if (rcd.rid > this->size()) return(std::vector<interval_type>());
		return(this->at(rcd.rid).findOverlapping(rcd.pos, rcd.pos + 1));
	}

private:
	void DedupeBlockList(void);

private:
	uint32_t n_intervals_;
    std::vector<std::string>   interval_strings_;
    std::vector< std::vector<interval_type> > interval_list_;
    std::vector<index_entry_type> block_list_;
    size_t  n_entries_; // equal to number of contigs
    pointer __entries;  // interval trees
};

}
}

#endif /* CONTAINERS_INTERVAL_CONTAINER_H_ */
