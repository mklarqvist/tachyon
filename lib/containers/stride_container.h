#ifndef CONTAINERS_STRIDE_CONTAINER_H_
#define CONTAINERS_STRIDE_CONTAINER_H_

#include <cassert>

#include "data_container.h"

namespace tachyon{
namespace containers{

#define YON_INTEGER_CONTAINER_DEFAULT_START_SIZE 1000

/**<
 * Primary container to handle integer data from data containers
 * This class should be considered for internal use only
 */
template <class return_primitive = U32>
class StrideContainer{
private:
	typedef StrideContainer   self_type;
    typedef std::size_t       size_type;
    typedef return_primitive  value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef DataContainer     data_container_type;

public:
    StrideContainer();
    StrideContainer(const size_type start_capacity);
    StrideContainer(const value_type uniform_value, const size_type n_entries);
    StrideContainer(const data_container_type& container);
    StrideContainer(const self_type& other);
    ~StrideContainer(void);

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

    /**<
     * Overloaded += operator for incrementally adding values
     * @param value Target value to be added
     */
    template <class T>
    inline void operator+=(const T& value){
    	if(this->size() + 1 == this->capacity())
    		this->resize();

    	this->__entries[this->n_entries] = value;
    }

    template <class T>
    inline void add(const T& value){ *this += value; }

    void resize(const size_type new_capacity){
    	if(new_capacity < this->capacity()){
    		this->n_entries = new_capacity;
    		return;
    	}

    	pointer old = this->__entries;
    	this->__entries = new value_type[new_capacity];
    	memcpy(this->data(), old, this->size()*sizeof(value_type));
    	delete [] old;
    	this->n_capacity =  new_capacity;
    }

    void resize(void){ this->resize(this->capacity()*2); }

private:
    /**<
     * Constructor invokes this function to in turn invoke
     * the correct `__allocate` funciton given the intrinsic
     * primitive
     * @param container Input data container
     */
    void __setup(const data_container_type& container);

    /**<
     * Called from `__allocate` to correctly copy data from
     * a given primitive to the current return primitive type
     * @param container Input data container
     */
    template <class intrinsic_type>
    void __allocate(const data_container_type& container);

    // Todo:
    bool determineUniformity(void);
    int findSmallestPrimitive(void);
    // 1) run find smallest primitive
    // 2) invoke stride container ctor with (larger_stride_container)

private:
    bool       isUniform_;
    size_type  n_capacity;
    size_type  n_entries;
    pointer    __entries;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer() :
	isUniform_(false),
	n_capacity(YON_INTEGER_CONTAINER_DEFAULT_START_SIZE),
	n_entries(0),
	__entries(new value_type[this->capacity()])
{
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const size_type start_capacity) :
	isUniform_(false),
	n_capacity(start_capacity),
	n_entries(0),
	__entries(new value_type[this->capacity()])
{
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const value_type uniform_value, const size_type n_entries) :
	isUniform_(true),
	n_capacity(n_entries),
	n_entries(n_entries),
	__entries(new value_type[this->capacity()])
{
	for(size_type i = 0; i < this->size(); ++i)
		this->__entries[i] = uniform_value;
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const data_container_type& container) :
	isUniform_(false),
	n_capacity(0),
	n_entries(0),
	__entries(nullptr)
{
	this->__setup(container);
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const self_type& other) :
	isUniform_(other.isUniform_),
	n_capacity(other.n_capacity),
	n_entries(other.n_entries),
	__entries(new value_type[this->size()])
{
	// Do not invoke memcpy as these two objects may have different primitive types
	for(size_type i = 0; i < this->size(); ++i)
		this->__entries[i] = other.__entries[i];
}

template <class return_primitive>
StrideContainer<return_primitive>::~StrideContainer(void){
	delete [] this->__entries;
}

template <class return_primitive>
void StrideContainer<return_primitive>::__setup(const data_container_type& container){
	switch(container.GetStridePrimitiveType()){
	case(YON_TYPE_8B):  this->__allocate<BYTE>(container); break;
	case(YON_TYPE_16B): this->__allocate<U16>(container);  break;
	case(YON_TYPE_32B): this->__allocate<U32>(container);  break;
	case(YON_TYPE_64B): this->__allocate<U64>(container);  break;
	case(YON_TYPE_FLOAT):
	case(YON_TYPE_DOUBLE):
	case(YON_TYPE_BOOLEAN):
	case(YON_TYPE_CHAR):
	case(YON_TYPE_STRUCT):
	case(YON_TYPE_UNKNOWN):
	default: std::cerr << utility::timestamp("ERROR") << "Illegal stride primitive: " << (int)container.header.stride_header.controller.type << std::endl; exit(1);
	}
}

template <class return_primitive>
template <class intrinsic_type>
void StrideContainer<return_primitive>::__allocate(const data_container_type& container){
	assert(container.buffer_strides_uncompressed.size() % sizeof(intrinsic_type) == 0);
	this->n_entries  = container.buffer_strides_uncompressed.size() / sizeof(intrinsic_type);
	this->__entries  = new value_type[this->size()];
	this->n_capacity = this->size();

	const intrinsic_type* const strides = reinterpret_cast<const intrinsic_type* const>(container.buffer_strides_uncompressed.data());

	for(size_type i = 0; i < this->size(); ++i)
		this->__entries[i] = strides[i];

	if(container.header.stride_header.controller.uniform)
		this->isUniform_ = true;
}

}
}



#endif /* CONTAINERS_STRIDE_CONTAINER_H_ */
