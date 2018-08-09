#ifndef CONTAINERS_CHECKSUM_CONTAINER_H_
#define CONTAINERS_CHECKSUM_CONTAINER_H_

#include "algorithm/digest/digest_manager.h"
#include "containers/variant_block.h"
#include "core/header/variant_header.h"
#include "containers/components/generic_iterator.h"

namespace tachyon{
namespace containers{

class ChecksumContainer {
private:
	typedef ChecksumContainer        self_type;
    typedef std::size_t              size_type;
    typedef algorithm::DigitalDigestPair  value_type;
    typedef value_type&              reference;
    typedef const value_type&        const_reference;
    typedef value_type*              pointer;
    typedef const value_type*        const_pointer;
    typedef io::BasicBuffer          buffer_type;
    typedef containers::VariantBlock block_type;
    typedef VariantHeader            header_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
	ChecksumContainer(void);
	ChecksumContainer(const size_type capacity);
	ChecksumContainer(const char* const data, const U32 l_data);
	ChecksumContainer(const buffer_type& buffer);
	~ChecksumContainer(void);

    // Element access
    inline reference at(const size_type& position){ return(this->__entries[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
    inline reference operator[](const size_type& position){ return(this->__entries[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
    inline pointer data(void){ return(this->__entries); }
    inline const_pointer data(void) const{ return(this->__entries); }
    inline reference front(void){ return(this->__entries[0]); }
    inline const_reference front(void) const{ return(this->__entries[0]); }
    inline reference back(void){ return(this->__entries[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries - 1]); }

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }
    inline const size_type& capacity(void) const{ return(this->n_capacity); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries]); }

    // Overload
    /**<
     * Allocates a number of empty elements
     * @param n_entries Number of entries to allocate
     * @return          Returns TRUE if successful or FALSE otherwise
     */
    bool allocate(const size_type n_entries);

    /**<
     * Add an element to the container
     * @param value Target element
     */
    void operator+=(const_reference value);

    /**<
     * Finalize all digital digest controllers in the container.
     * This functions has to be called prior to reading/writing to disk
     */
    void finalize(void);

    /**<
     * Updates all target digital digest records in the container given
     * the input data block
     * @param block    Input data block container
     * @param header   Header
     * @return         Returns TRUE upon success or FALSE otherwise
     */
	bool update(const block_type& block, const header_type& header);

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& container){
		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries), sizeof(size_type));
		for(size_type i = 0; i < container.size(); ++i)
			out << container[i];

		return(out);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container){
		stream.read((char*)reinterpret_cast<size_type*>(&container.n_entries), sizeof(size_type));
		container.allocate(container.n_entries);
		for(size_type i = 0; i < container.size(); ++i)
			stream >> container[i];

		return(stream);
	}

private:
    size_t  n_entries;
    size_t  n_capacity;
    pointer __entries;
};

}
}


#endif /* CONTAINERS_CHECKSUM_CONTAINER_H_ */
