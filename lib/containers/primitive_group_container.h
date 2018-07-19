#ifndef CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_
#define CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_

#include "primitive_container.h"

namespace tachyon{
namespace containers{

template <class return_type>
class PrimitiveGroupContainer{
private:
    typedef PrimitiveGroupContainer self_type;
    typedef PrimitiveContainer<return_type> value_type;
    typedef std::size_t             size_type;
    typedef value_type&             reference;
    typedef const value_type&       const_reference;
    typedef value_type*             pointer;
    typedef const value_type*       const_pointer;
    typedef std::ptrdiff_t          difference_type;
    typedef DataContainer           data_container_type;

public:
    PrimitiveGroupContainer();
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
	inline bool empty(void) const{ return(this->__n_objects == 0); }
	inline const size_type& size(void) const{ return(this->__n_objects); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__containers[0]); }
	inline iterator end(){ return iterator(&this->__containers[this->__n_objects]); }
	inline const_iterator begin() const{ return const_iterator(&this->__containers[0]); }
	inline const_iterator end() const{ return const_iterator(&this->__containers[this->__n_objects]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->__containers[this->__n_objects]); }

private:
	template <class actual_primitive_type>
	void __setup(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each);

private:
    size_type __n_objects;
    pointer   __containers;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer() : __n_objects(0), __containers(nullptr){}

template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each) :
	__n_objects(n_objects),
	__containers(static_cast<pointer>(::operator new[](this->__n_objects*sizeof(value_type))))
{

	if(container.header.data_header.isSigned()){
		switch(container.header.data_header.getPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setup<SBYTE>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_16B):    (this->__setup<S16>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_32B):    (this->__setup<S32>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_64B):    (this->__setup<S64>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << "Disallowed: " << container.header.data_header.getPrimitiveType() << std::endl; return;
		}
	} else {
		switch(container.header.data_header.getPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setup<BYTE>(container, offset, n_objects, strides_each));   break;
		case(YON_TYPE_16B):    (this->__setup<U16>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_32B):    (this->__setup<U32>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_64B):    (this->__setup<U64>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << "Disallowed: " << container.header.data_header.getPrimitiveType() << std::endl; return;
		}
	}
}

template <class return_type>
PrimitiveGroupContainer<return_type>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->__n_objects; ++i)
		((this->__containers + i)->~value_type)();

	::operator delete[](static_cast<void*>(this->__containers));
}

template <class return_type>
template <class actual_primitive_type>
void PrimitiveGroupContainer<return_type>::__setup(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each){
	U32 current_offset = offset;
	for(U32 i = 0; i < this->__n_objects; ++i){
		new( &this->__containers[i] ) value_type( container, current_offset, strides_each );
		current_offset += strides_each * sizeof(actual_primitive_type);
	}
}

}
}



#endif /* CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_ */
