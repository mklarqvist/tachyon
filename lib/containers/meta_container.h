#ifndef CONTAINERS_META_CONTAINER_H_
#define CONTAINERS_META_CONTAINER_H_

#include "variant_block.h"
#include "core/meta_entry.h"
#include "components/generic_iterator.h"

namespace tachyon{
namespace containers{

class MetaContainer {
public:
	typedef MetaContainer      self_type;
    typedef std::size_t        size_type;
    typedef core::MetaEntry    value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef VariantBlock       block_type;
    typedef VariantBlockHeader block_header_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
	MetaContainer(const block_type& block);
	MetaContainer(const self_type& other);
	MetaContainer(self_type&& other) noexcept;
	MetaContainer& operator=(const self_type& other);
	MetaContainer& operator=(self_type&& other) noexcept;
	~MetaContainer(void);

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

    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries]); }

    /**<
     * Convert relative offsets into absolute offsets.
     * @param header Reference block header providing the relative minimum position.
     */
    void CorrectRelativePositions(const block_header_type& header){
    	for(uint32_t i = 0; i < this->size(); ++i)
    		this->at(i).position += header.minPosition;
    }

private:
    void Setup(const block_type& block);

private:
    size_t  n_entries;
    pointer __entries;
};

}
}

#endif /* CONTAINERS_META_CONTAINER_H_ */
