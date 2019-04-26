#ifndef ALGORITHM_DigitalDigest_H_
#define ALGORITHM_DigitalDigest_H_

#include <cstring>
#include <openssl/sha.h>

#include "buffer.h"
#include "digest.h"

namespace tachyon{
namespace algorithm{

class DigestManager{
private:
	typedef DigestManager     self_type;

protected:
	typedef DigitalDigestPair value_type;
	typedef DigitalDigest     digest_type;
	typedef std::size_t       size_type;
	typedef value_type&       reference;
	typedef const value_type& const_reference;
	typedef value_type*       pointer;
	typedef const value_type* const_pointer;

public:
	DigestManager() :
		n_entries_(0),
		n_capacity_(100),
		__entries(new value_type[this->n_capacity_])
	{}

	DigestManager(const size_type start_capacity) :
		n_entries_(start_capacity),
		n_capacity_(start_capacity),
		__entries(new value_type[this->n_capacity_])
	{}

	DigestManager(const self_type& other) :
		n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		__entries(new value_type[this->n_capacity_])
	{
		for (uint32_t i = 0; i < this->size(); ++i) this->__entries[i] = other.__entries[i];
	}

	virtual ~DigestManager() { delete [] this->__entries; }

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const { return *ptr_; }
		pointer operator->() const { return ptr_; }
		bool operator==(const self_type& rhs) const { return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const { return ptr_ != rhs.ptr_; }
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
		const_reference operator*() const { return *ptr_; }
		const_pointer operator->() const { return ptr_; }
		bool operator==(const self_type& rhs) const { return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const { return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position) { return(this->__entries[position]); }
	inline const_reference at(const size_type& position) const { return(this->__entries[position]); }
	inline reference operator[](const size_type& position) { return(this->__entries[position]); }
	inline const_reference operator[](const size_type& position) const { return(this->__entries[position]); }
	inline pointer data(void) { return(this->__entries); }
	inline const_pointer data(void) const { return(this->__entries); }
	inline reference front(void) { return(this->__entries[0]); }
	inline const_reference front(void) const { return(this->__entries[0]); }
	inline reference back(void) { return(this->__entries[this->n_entries_ - 1]); }
	inline const_reference back(void) const { return(this->__entries[this->n_entries_ - 1]); }

	// Capacity
	inline bool empty(void) const { return(this->n_entries_ == 0); }
	inline const size_type& size(void) const { return(this->n_entries_); }
	inline const size_type& capacity(void) const { return(this->n_capacity_); }

	// Iterator
	inline iterator begin() { return iterator(&this->__entries[0]); }
	inline iterator end() { return iterator(&this->__entries[this->n_entries_]); }
	inline const_iterator begin() const { return const_iterator(&this->__entries[0]); }
	inline const_iterator end() const { return const_iterator(&this->__entries[this->n_entries_]); }
	inline const_iterator cbegin() const { return const_iterator(&this->__entries[0]); }
	inline const_iterator cend() const { return const_iterator(&this->__entries[this->n_entries_]); }

	void finalize(void) {
		for (uint32_t i = 0; i < this->size(); ++i) this->at(i).finalize();
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& container) {
		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_), sizeof(size_type));
		for (size_type i = 0; i < container.size(); ++i) out << container[i];

		return(out);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container) {
		stream.read((char*)reinterpret_cast<size_type*>(&container.n_entries_), sizeof(size_type));
		delete [] container.__entries;
		container.__entries = new value_type[container.n_entries_];
		for (size_type i = 0; i < container.size(); ++i) stream >> container[i];

		return(stream);
	}

protected:
	size_type n_entries_;
	size_type n_capacity_;
	pointer   __entries;
};

}
}

#endif /* ALGORITHM_DigitalDigest_H_ */
