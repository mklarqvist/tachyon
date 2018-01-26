#ifndef CONTAINERS_PRIMITIVEGROUPCONTAINER_H_
#define CONTAINERS_PRIMITIVEGROUPCONTAINER_H_

namespace tachyon{
namespace core{

template <class T>
class PrimitiveGroupContainer{
private:
    typedef PrimitiveGroupContainer self_type;
    typedef PrimitiveContainer<T> value_type;
    typedef std::size_t           size_type;
    typedef value_type&           reference;
    typedef const value_type&     const_reference;
    typedef value_type*           pointer;
    typedef const value_type*     const_pointer;
    typedef std::ptrdiff_t        difference_type;
    typedef DataContainer         data_container_type;

public:
    // PrimitiveGroupContainer();
    PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each);
    ~PrimitiveGroupContainer(void);

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
	inline reference at(const size_type& position){ return(this->__containers[position]); }
	inline const_reference at(const size_type& position) const{ return(this->__containers[position]); }
	inline reference operator[](const size_type& position){ return(this->__containers[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->__containers[position]); }
	inline pointer data(void){ return(this->__containers); }
	inline const_pointer data(void) const{ return(this->__containers); }
	inline reference front(void){ return(this->__containers[0]); }
	inline const_reference front(void) const{ return(this->__containers[0]); }
	inline reference back(void){ return(this->__containers[this->__n_objects - 1]); }
	inline const_reference back(void) const{ return(this->__containers[this->__n_objects - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->__n_objects == 0); }
	inline const size_type& size(void) const{ return(this->__n_objects); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__containers[0]); }
	inline iterator end(){ return iterator(&this->__containers[this->__n_objects - 1]); }
	inline const_iterator begin() const{ return const_iterator(&this->__containers[0]); }
	inline const_iterator end() const{ return const_iterator(&this->__containers[this->__n_objects - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->__containers[this->__n_objects - 1]); }

private:
    size_type __n_objects;
    pointer   __containers;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class T>
PrimitiveGroupContainer<T>::PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each) :
	__n_objects(n_objects),
	__containers(static_cast<pointer>(::operator new[](this->__n_objects*sizeof(value_type))))
{
	U32 current_offset = offset;
	for(U32 i = 0; i < this->__n_objects; ++i){
		new( &this->__containers[i] ) value_type( container, current_offset, strides_each );
		current_offset += strides_each * sizeof(T);
	}
}

template <class T>
PrimitiveGroupContainer<T>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->__n_objects; ++i)
		((this->__containers + i)->~value_type)();

	::operator delete[](static_cast<void*>(this->__containers));
}

}
}



#endif /* CONTAINERS_PRIMITIVEGROUPCONTAINER_H_ */
