#ifndef CONTAINERS_STRIDE_CONTAINER_H_
#define CONTAINERS_STRIDE_CONTAINER_H_

#include <cassert>

#include "components/generic_iterator.h"
#include "data_container.h"

namespace tachyon{
namespace containers{

/**<
 * Do not use resize or operator+= functionality external to this
 * class as the StrideContainer is a contextual representation of
 * the strides stored in the VariantBlock structure. As such, any
 * changes to the strides should be represented there, not in
 * this representation.
 */
template <class return_primitive = uint32_t>
class StrideContainer {
public:
	typedef StrideContainer   self_type;
    typedef std::size_t       size_type;
    typedef return_primitive  value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef DataContainer     data_container_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    StrideContainer();
    StrideContainer(const size_type start_capacity);
    StrideContainer(const value_type uniform_value, const size_type n_entries_);
    StrideContainer(const data_container_type& container);
    StrideContainer(const self_type& other);
    ~StrideContainer(void);
    StrideContainer(self_type&& other) noexcept;
    StrideContainer& operator=(const self_type& other);
    StrideContainer& operator=(self_type&& other) noexcept;

    // Element access
    inline reference at(const size_type& position){ return(this->entries_[position]); }
    inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
    inline reference operator[](const size_type& position){ return(this->entries_[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
    inline pointer data(void){ return(this->entries_); }
    inline const_pointer data(void) const{ return(this->entries_); }
    inline reference front(void){ return(this->entries_[0]); }
    inline const_reference front(void) const{ return(this->entries_[0]); }
    inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

    // Capacity
    inline bool& IsUniform(void) const{ return(this->is_uniform_); }
    inline bool empty(void) const{ return(this->n_entries_ == 0); }
    inline const size_type& size(void) const{ return(this->n_entries_); }
    inline const size_type& capacity(void) const{ return(this->n_capacity_); }

    // Iterator
    inline iterator begin(){ return iterator(&this->entries_[0]); }
    inline iterator end(){ return iterator(&this->entries_[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries_]); }

    /**<
     * Overloaded += operator for incrementally adding values to
     * this container. External uses of this function should be
     * avoided as StrideContainer objects are not used when writing
     * objects to disk. Data written from disk is sourced from the
     * DataBlock buffers themselves.
     * @param value Target value to be added
     */
    template <class T>
    inline void operator+=(const T& value){
    	if(this->size() + 1 == this->capacity())
    		this->resize();

    	this->entries_[this->n_entries_] = value;
    }

    template <class T>
    inline void add(const T& value){ *this += value; }

    void resize(const size_type new_capacity);
    inline void resize(void){ this->resize(this->capacity()*2); }

private:
    /**<
     * Constructor invokes this function to in turn invoke
     * the correct Allocate() function given the intrinsic
     * primitive.
     * @param container Input data container
     */
    void Setup(const data_container_type& container);

    /**<
     * Called from Setup() to correctly copy data from
     * a given primitive to the desired return primitive type.
     * Note that this function is double templated: the first
     * correspond to the desired return primitive type and the
     * second the primitive type used to store the data in the
     * source byte stream.
     * @param container Input data container.
     */
    template <class intrinsic_type>
    void Allocate(const data_container_type& container);


private:
    bool       is_uniform_;
    size_type  n_capacity_;
    size_type  n_entries_;
    pointer    entries_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer() :
	is_uniform_(false),
	n_capacity_(1000),
	n_entries_(0),
	entries_(new value_type[this->capacity()])
{
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const size_type start_capacity) :
	is_uniform_(false),
	n_capacity_(start_capacity),
	n_entries_(0),
	entries_(new value_type[this->capacity()])
{
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const value_type uniform_value, const size_type n_entries_) :
	is_uniform_(true),
	n_capacity_(n_entries_),
	n_entries_(n_entries_),
	entries_(new value_type[this->capacity()])
{
	for(size_type i = 0; i < this->size(); ++i)
		this->entries_[i] = uniform_value;
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const data_container_type& container) :
	is_uniform_(false),
	n_capacity_(0),
	n_entries_(0),
	entries_(nullptr)
{
	this->Setup(container);
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(const self_type& other) :
	is_uniform_(other.is_uniform_),
	n_capacity_(other.n_capacity_),
	n_entries_(other.n_entries_),
	entries_(new value_type[this->size()])
{
	// Do not use a memcpy call here as the stored primitive type
	// could be different from the desired return primitive type.
	for(size_type i = 0; i < this->size(); ++i)
		this->entries_[i] = other.entries_[i];
}

template <class return_primitive>
StrideContainer<return_primitive>::StrideContainer(self_type&& other) noexcept :
	is_uniform_(other.is_uniform_),
	n_capacity_(other.n_capacity_),
	n_entries_(other.n_entries_),
	entries_(nullptr)
{
	std::swap(this->entries_, other.entries_);
}

template <class return_primitive>
StrideContainer<return_primitive>::~StrideContainer(void){
	delete [] this->entries_;
}

template <class return_primitive>
StrideContainer<return_primitive>& StrideContainer<return_primitive>::operator=(self_type&& other) noexcept {
	this->is_uniform_ = other.is_uniform_;
	this->n_capacity_ = other.n_capacity_;
	this->n_entries_  = other.n_entries_;
	delete [] this->entries_; this->entries_ = nullptr;
	std::swap(this->entries_, other.entries_);
}

template <class return_primitive>
StrideContainer<return_primitive>& StrideContainer<return_primitive>::operator=(const self_type& other) {
	this->is_uniform_ = other.is_uniform_;
	this->n_capacity_ = other.n_capacity_;
	this->n_entries_  = other.n_entries_;
	delete [] this->entries_;
	this->entries_ = new value_type[this->capacity()];
	memcpy(this->entries_, other.entries_, sizeof(value_type)*this->size());
}

template <class return_primitive>
void StrideContainer<return_primitive>::Setup(const data_container_type& container){
	switch(container.GetStridePrimitiveType()){
	case(YON_TYPE_8B):  this->Allocate<uint8_t>(container);  break;
	case(YON_TYPE_16B): this->Allocate<uint16_t>(container); break;
	case(YON_TYPE_32B): this->Allocate<uint32_t>(container); break;
	case(YON_TYPE_64B): this->Allocate<uint64_t>(container); break;
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
void StrideContainer<return_primitive>::Allocate(const data_container_type& container){
	// Assert that the input type is divisible by the primitive
	// type byte width. If this is not true then the data it
	// guaranteed to be corrupted. Alternatively, the primitive
	// type could be misrepresented in the data header.
	assert(container.strides_uncompressed.size() % sizeof(intrinsic_type) == 0);

	this->n_entries_  = container.strides_uncompressed.size() / sizeof(intrinsic_type);
	this->entries_    = new value_type[this->size()];
	this->n_capacity_ = this->size();

	// Cast the buffer as the primitive type it was stored as. This
	// is required for correct interpretation of the byte stream.
	const intrinsic_type* const strides = reinterpret_cast<const intrinsic_type* const>(container.strides_uncompressed.data());

	// Iterate over available data and copy it over. Do not use a memcpy
	// call here as the stored primitive type could be different from the
	// desired return primitive type.
	for(size_type i = 0; i < this->size(); ++i)
		this->entries_[i] = strides[i];

	// If the data is uniform then trigger a flag remembering
	// this fact.
	if(container.header.stride_header.controller.uniform)
		this->is_uniform_ = true;
}

template <class return_primitive>
void StrideContainer<return_primitive>::resize(const size_type new_capacity){
	if(new_capacity < this->capacity()){
		this->n_entries_ = new_capacity;
		return;
	}

	pointer old    = this->entries_;
	this->entries_ = new value_type[new_capacity];
	memcpy(this->data(), old, this->size()*sizeof(value_type));
	delete [] old;
	this->n_capacity_ =  new_capacity;
}

}
}



#endif /* CONTAINERS_STRIDE_CONTAINER_H_ */
