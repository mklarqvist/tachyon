#ifndef CONTAINERS_ABSTRACTINTEGERCONTAINER_H_
#define CONTAINERS_ABSTRACTINTEGERCONTAINER_H_

#include "../iterator/IteratorIntegerReference.h"

namespace Tachyon{
namespace Core{

template <class return_primitive>
class AbstractIntegerContainer{
    typedef Iterator::IteratorIntegerReference<return_primitive> value_type;
    typedef std::size_t       size_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;

public:
    AbstractIntegerContainer(const Container& container);
    ~AbstractIntegerContainer();

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
	inline reference at(const size_type& position){ return(this->__iterators[position]); }
	inline const_reference at(const size_type& position) const{ return(this->__iterators[position]); }
	inline reference operator[](const size_type& position){ return(this->__iterators[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->__iterators[position]); }
	inline pointer data(void){ return(this->__iterators); }
	inline const_pointer data(void) const{ return(this->__iterators); }
	inline reference front(void){ return(this->__iterators[0]); }
	inline const_reference front(void) const{ return(this->__iterators[0]); }
	inline reference back(void){ return(this->__iterators[this->n_entries - 1]); }
	inline const_reference back(void) const{ return(this->__iterators[this->n_entries - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__iterators[0]); }
	inline iterator end(){ return iterator(&this->__iterators[this->n_entries - 1]); }
	inline const_iterator begin() const{ return const_iterator(&this->__iterators[0]); }
	inline const_iterator end() const{ return const_iterator(&this->__iterators[this->n_entries - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__iterators[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->__iterators[this->n_entries - 1]); }

private:
    size_t      n_entries;   // number of iterators = number of stride primitives
    char*       __buffer;    // copy data
    value_type* __iterators; // iterators
};

template <class return_primitive>
AbstractIntegerContainer<return_primitive>::AbstractIntegerContainer(const Container& container) :
	n_entries(0),
	__buffer(new char[container.buffer_data_uncompressed.size()]),
	__iterators(nullptr)
{
	if(container.buffer_data_uncompressed.size() == 0)
		return;

	switch(container.header_stride.controller.type){
	case(Core::YON_TYPE_8B):  this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(BYTE); break;
	case(Core::YON_TYPE_16B): this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U16);  break;
	case(Core::YON_TYPE_32B): this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U32);  break;
	case(Core::YON_TYPE_64B): this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U64);  break;
	}

	if(this->n_entries == 0)
		return;

	//
	this->__iterators = new value_type[this->n_entries];
}

}
}



#endif /* CONTAINERS_ABSTRACTINTEGERCONTAINER_H_ */
