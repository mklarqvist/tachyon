#ifndef BCF_BCFREADER_H_
#define BCF_BCFREADER_H_

#include <cassert>

#include "BCFEntry.h"
#include "../compression/BGZFController.h"


namespace tachyon {
namespace bcf {

class BCFReader{
    typedef BCFReader          self_type;
    typedef BCFEntry           value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef std::ptrdiff_t     difference_type;
    typedef std::size_t        size_type;
    typedef io::BasicBuffer    buffer_type;
    typedef io::BGZFController bgzf_controller_type;
    typedef vcf::VCFHeader     header_type;
	typedef core::HeaderContig contig_type;

public:
	enum bcf_reader_state{BCF_INIT, BCF_OK, BCF_ERROR, BCF_EOF, BCF_STREAM_ERROR};

public:
	BCFReader();
	BCFReader(const std::string& file_name);
	~BCFReader();

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
	inline reference at(const size_type& position){ return(this->entries[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries[position]); }
	inline reference operator[](const size_type& position){ return(this->entries[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries[position]); }
	inline pointer data(void){ return(this->entries); }
	inline const_pointer data(void) const{ return(this->entries); }
	inline reference front(void){ return(this->entries[0]); }
	inline const_reference front(void) const{ return(this->entries[0]); }
	inline reference back(void){ return(this->entries[this->n_entries - 1]); }
	inline const_reference back(void) const{ return(this->entries[this->n_entries - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline const size_type& capacity(void) const{ return(this->n_capacity); }

	// Iterator
	inline iterator begin(){ return iterator(&this->entries[0]); }
	inline iterator end()  { return iterator(&this->entries[this->n_entries - 1]); }
	inline const_iterator begin()  const{ return const_iterator(&this->entries[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->entries[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->entries[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->entries[this->n_entries - 1]); }

	/**<
	 * Attempts to open a target input file. Internally
	 * checks if the input file is an actual BCF file and
	 * the first TGZF can be opened and the BCF header is
	 * valid.
	 * @param input Input target BCF file
	 * @return      Returns TRUE upon success or FALSE otherwise
	 */
	bool open(const std::string input);
	bool open(void);

	/**<
	 * Loads another TGZF block into memory
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool nextBlock(void);

	/**<
	 * Attempts to overload `entry` input BCFEntry with
	 * data
	 * @param entry Input BCFEntry that will be overloaded
	 * @return      Returns TRUE upon success or FALSE otherwise
	 */
	bool nextVariant(reference entry);

	/**<
	 * Attempts to load either `n_variants` number of variants or
	 * variants covering >= `bp_window` base pairs. The function
	 * is successful whenever n_variants or bp_window is satisfied
	 * or if there is no more variants to load.
	 * @param n_variants     Number of variants
	 * @param bp_window      Non-overlapping window size in base-pairs
	 * @param across_contigs Allow the algorithm to span over two or more different chromosomes
	 * @return               Returns TRUE upon success or FALSE otherwise
	 */
	bool getVariants(const U32 n_variants, const double bp_window, bool across_contigs = false); // get N number of variants into buffer

	inline const bool hasCarryOver(void) const{ return(this->n_carry_over); }

private:
	/**<
	 * Parse the TGZF header of a block given
	 * the current buffer data loaded.
	 * Internal use only
	 * @return Returns TRUE on success or FALSE otherwise
	 */
	bool parseHeader(void);

public:
	std::string          file_name;
	std::ifstream        stream;
	U64                  filesize;
	U32                  current_pointer;
	buffer_type          buffer;
	buffer_type          header_buffer;
	bgzf_controller_type bgzf_controller;
	header_type          header;
	bcf_reader_state     state;
	size_type            n_entries;
	size_type            n_capacity;
	size_type            n_carry_over;
	pointer              entries;
};

}
}

#endif /* BCF_BCFREADER_H_ */
