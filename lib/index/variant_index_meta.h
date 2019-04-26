#ifndef INDEX_INDEX_META_CONTAINER_H_
#define INDEX_INDEX_META_CONTAINER_H_

#include "variant_index_meta_entry.h"
#include "generic_iterator.h"

namespace tachyon{
namespace index{

class VariantIndexMeta{
public:
	typedef VariantIndexMeta   self_type;
    typedef std::size_t        size_type;
    typedef VariantIndexMetaEntry value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    VariantIndexMeta() :
		n_entries(0),
		n_capacity(1000),
		__entries(new value_type[this->n_capacity])
	{}

    VariantIndexMeta(const self_type& other) :
    	n_entries(other.n_entries),
		n_capacity(other.n_capacity),
		__entries(new value_type[this->n_capacity])
    {
    	for (uint32_t i = 0; i < this->size(); ++i)
    		this->__entries[i] = other.__entries[i];
    }

	~VariantIndexMeta() {
		delete [] this->__entries;
	}

	// Element access
	inline reference at(const size_type& position) { return(this->__entries[position]); }
	inline const_reference at(const size_type& position) const { return(this->__entries[position]); }
	inline reference operator[](const size_type& position) { return(this->__entries[position]); }
	inline const_reference operator[](const size_type& position) const { return(this->__entries[position]); }
	inline pointer data(void) { return(this->__entries); }
	inline const_pointer data(void) const { return(this->__entries); }
	inline reference front(void) { return(this->__entries[0]); }
	inline const_reference front(void) const { return(this->__entries[0]); }
	inline reference back(void) { return(this->__entries[this->n_entries - 1]); }
	inline const_reference back(void) const { return(this->__entries[this->n_entries - 1]); }

	// Capacity
	inline bool empty(void) const { return(this->n_entries == 0); }
	inline const size_type& size(void) const { return(this->n_entries); }
	inline const size_type& capacity(void) const { return(this->n_capacity); }

	// Iterator
	inline iterator begin() { return iterator(&this->__entries[0]); }
	inline iterator end() { return iterator(&this->__entries[this->n_entries]); }
	inline const_iterator begin() const { return const_iterator(&this->__entries[0]); }
	inline const_iterator end() const { return const_iterator(&this->__entries[this->n_entries]); }
	inline const_iterator cbegin() const { return const_iterator(&this->__entries[0]); }
	inline const_iterator cend() const { return const_iterator(&this->__entries[this->n_entries]); }

	inline self_type& operator+=(const_reference index_entry) {
		if (this->size() + 1 == this->n_capacity)
			this->resize();


		this->__entries[this->n_entries++] = index_entry;
		return(*this);
	}
	inline self_type& add(const_reference index_entry) { return(*this += index_entry); }

	void resize(void) {
		pointer temp = this->__entries;

		this->n_capacity *= 2;
		this->__entries = new value_type[this->capacity()];

		// Lift over values from old addresses
		for (uint32_t i = 0; i < this->size(); ++i)
			this->__entries[i] = temp[i];

		delete [] temp;
	}

	void resize(const uint32_t new_size) {
		pointer temp = this->__entries;

		this->n_capacity = new_size;
		this->__entries = new value_type[this->capacity()];

		// Lift over values from old addresses
		for (uint32_t i = 0; i < this->size(); ++i)
			this->__entries[i] = temp[i];

		delete [] temp;
	}

	void reserve(const uint32_t new_size) {
		this->resize(new_size);
		this->n_entries = new_size;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry) {
		stream.write(reinterpret_cast<const char*>(&entry.n_entries), sizeof(size_type));
		for (uint32_t i = 0; i < entry.size(); ++i) stream << entry[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry) {
		stream.read(reinterpret_cast<char*>(&entry.n_entries), sizeof(size_type));

		if (entry.size() > entry.capacity()) {
			delete [] entry.__entries;
			entry.n_capacity = entry.size() + 100;
			entry.__entries = new value_type[entry.capacity()];
		}

		for (uint32_t i = 0; i < entry.size(); ++i)
			stream >> entry[i];

		return(stream);
	}

private:
	size_type n_entries;
	size_type n_capacity;
	pointer   __entries;
};


}
}



#endif /* INDEX_INDEX_META_CONTAINER_H_ */
