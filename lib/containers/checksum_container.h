#ifndef CONTAINERS_CHECKSUM_CONTAINER_H_
#define CONTAINERS_CHECKSUM_CONTAINER_H_

#include "algorithm/digest/digest_manager.h"
#include "header_footer.h"
#include "generic_iterator.h"
#include "variant_container.h"

namespace tachyon{
namespace containers{

class ChecksumContainer {
public:
	typedef ChecksumContainer        self_type;
    typedef std::size_t              size_type;
    typedef algorithm::DigitalDigestPair  value_type;
    typedef value_type&              reference;
    typedef const value_type&        const_reference;
    typedef value_type*              pointer;
    typedef const value_type*        const_pointer;
    typedef yon_buffer_t             buffer_type;
    typedef yon1_vb_t                block_type;
    typedef yon_vnt_hdr_t            header_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
	ChecksumContainer(void);
	ChecksumContainer(const size_type capacity);
	ChecksumContainer(const char* const data, const uint32_t l_data);
	ChecksumContainer(const buffer_type& buffer);
	~ChecksumContainer(void);

    // Element access
    inline reference at(const size_type& position) { return(this->entries_[position]); }
    inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
    inline reference operator[](const size_type& position) { return(this->entries_[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
    inline pointer data(void) { return(this->entries_); }
    inline const_pointer data(void) const{ return(this->entries_); }
    inline reference front(void) { return(this->entries_[0]); }
    inline const_reference front(void) const{ return(this->entries_[0]); }
    inline reference back(void) { return(this->entries_[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->entries_[this->n_entries - 1]); }

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const size_type& capacity(void) const{ return(this->n_capacity); }

    // Iterator
    inline iterator begin() { return iterator(&this->entries_[0]); }
    inline iterator end() { return iterator(&this->entries_[this->n_entries]); }
    inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries]); }

    // Overload
    /**<
     * Allocates a number of empty elements
     * @param n_entries Number of entries to allocate
     * @return          Returns TRUE if successful or FALSE otherwise
     */
    bool Allocate(const size_type n_entries);

    /**<
     * Add an element to the container
     * @param value Target element
     */
    void operator+=(const_reference value);

    /**<
     * Finalize all digital digest controllers in the container.
     * This functions has to be called prior to reading/writing to disk
     */
    void Finalize(void);

    /**<
     * Updates all target digital digest records in the container given
     * the input data block
     * @param block    Input data block container
     * @param header   Header
     * @return         Returns TRUE upon success or FALSE otherwise
     */
	bool Update(const block_type& block, const header_type& header);

	friend std::ostream& operator<<(std::ostream& out, const self_type& container);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container);

private:
    size_t  n_entries;
    size_t  n_capacity;
    pointer entries_;
};

}
}


#endif /* CONTAINERS_CHECKSUM_CONTAINER_H_ */
