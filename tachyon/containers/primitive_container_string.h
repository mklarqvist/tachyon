#ifndef CONTAINERS_PRIMITIVE_CONTAINER_STRING_H_
#define CONTAINERS_PRIMITIVE_CONTAINER_STRING_H_

#include "primitive_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

/**<
 * Primitive container for strings
 */
template <>
class PrimitiveContainer<std::string>{
private:
    typedef std::size_t          size_type;
    typedef std::string          value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef StrideContainer<U32> stride_container_type;
    typedef DataContainer        data_container_type;

public:
    PrimitiveContainer();
    PrimitiveContainer(const data_container_type& container);
    ~PrimitiveContainer(void);

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
    inline reference back(void){ return(this->__entries[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries - 1]); }

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries]); }

private:
    void __setup(const data_container_type& container);
    void __setup(const data_container_type& container, const U32 stride_size);

private:
    size_t  n_entries;
    pointer __entries;
};


// IMPLEMENTATION -------------------------------------------------------------


PrimitiveContainer<std::string>::PrimitiveContainer() :
	n_entries(0),
	__entries(nullptr)
{

}

PrimitiveContainer<std::string>::PrimitiveContainer(const data_container_type& container) :
	n_entries(n_entries),
	__entries(new value_type[n_entries])
{
	if(container.header.getPrimitiveType() != tachyon::core::YON_TYPE_CHAR){
		std::cerr << "disallowed type" << std::endl;
		return;
	}

	if(container.header.hasMixedStride()){
		this->__setup(container);
	} else {
		this->__setup(container, container.header.stride);
	}
}

PrimitiveContainer<std::string>::~PrimitiveContainer(void){
	delete [] this->__entries;
}

void PrimitiveContainer<std::string>::__setup(const data_container_type& container){
	//const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.buffer_data_uncompressed.buffer[offset]);
	//for(U32 i = 0; i < this->n_entries; ++i)
	//	this->__entries[i] = data[i];
}

void PrimitiveContainer<std::string>::__setup(const data_container_type& container, const U32 stride_size){
	//const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.buffer_data_uncompressed.buffer[offset]);
	//for(U32 i = 0; i < this->n_entries; ++i)
	//	this->__entries[i] = data[i];
}

}
}



#endif /* CONTAINERS_PRIMITIVE_CONTAINER_STRING_H_ */
