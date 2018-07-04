#ifndef CONTAINERS_INTERVAL_CONTAINER_H_
#define CONTAINERS_INTERVAL_CONTAINER_H_

#include <regex>

#include "support/magic_constants.h"
#include "support/type_definitions.h"
#include "support/helpers.h"
#include "third_party/intervalTree.h"
#include "index/index.h"
#include "core/header/variant_header.h"
#include "core/meta_entry.h"

namespace tachyon{
namespace containers{

class IntervalContainer {
private:
	typedef IntervalContainer      self_type;
    typedef std::size_t            size_type;
    typedef algorithm::Interval<U32, S64> interval_type;
    typedef algorithm::IntervalTree<U32, S64>  value_type;
    typedef value_type&            reference;
    typedef const value_type&      const_reference;
    typedef value_type*            pointer;
    typedef const value_type*      const_pointer;
    typedef index::Index           index_type;
    typedef index::IndexEntry      index_entry_type;
    typedef core::VariantHeader    header_type;
    typedef core::MetaEntry        meta_entry_type;

public:
    IntervalContainer();
    ~IntervalContainer(void);

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

    // Element access
    inline reference at(const size_type& position){ return(this->__entries[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
    inline reference operator[](const size_type& position){ return(this->__entries[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
    inline pointer data(void){ return(this->__entries); }
    inline const_pointer data(void) const{ return(this->__entries); }
    inline reference front(void){ return(this->__entries[0]); }
    inline const_reference front(void) const{ return(this->__entries[0]); }
    inline reference back(void){ return(this->__entries[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries_ - 1]); }

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries_ == 0); }
    inline const size_type& size(void) const{ return(this->n_entries_); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries_]); }

	inline const bool hasBlockList(void) const{ return(this->block_list_.size()); }
	inline std::vector<index_entry_type>& getBlockList(void){ return(this->block_list_); }
	inline const std::vector<index_entry_type>& getBlockList(void) const{ return(this->block_list_); }
	inline const bool hasIntervals(void) const{ return(this->interval_list_.size()); }

    // Interpret
    bool validateIntervalStrings(std::vector<std::string>& interval_strings);

    /**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool parseIntervals(std::vector<std::string>& interval_strings, const header_type& header, const index_type& index);

	bool build(const header_type& header);

	inline std::vector<interval_type> find_overlaps(const U32& contigID, const S64& start_position, const S64& end_position) const{
		if(contigID > this->size()) return(std::vector<interval_type>());
		return(this->at(contigID).findOverlapping(start_position, end_position));
	}

	inline std::vector<interval_type> find_overlaps(const meta_entry_type& meta_entry) const{
		if(meta_entry.contigID > this->size()) return(std::vector<interval_type>());
		return(this->at(meta_entry.contigID).findOverlapping(meta_entry.position, meta_entry.position + 1));
	}

private:
	void dedupeBlockList(void);

private:
	U32 n_intervals_;
    std::vector<std::string>   interval_strings_;
    std::vector< std::vector<interval_type> > interval_list_;
    std::vector<index_entry_type> block_list_;
    size_t  n_entries_; // equal to number of contigs
    pointer __entries;  // interval trees
};

}
}

#endif /* CONTAINERS_INTERVAL_CONTAINER_H_ */
