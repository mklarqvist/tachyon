#ifndef CONTAINER_INFOCONTAINER_H_
#define CONTAINER_INFOCONTAINER_H_

#include "datacontainer.h"
#include "primitive_container.h"

namespace tachyon{
namespace containers{

template <class return_type>
class InfoContainer{
private:
    typedef InfoContainer                   self_type;
    typedef PrimitiveContainer<return_type> value_type;
    typedef value_type&                     reference;
    typedef const value_type&               const_reference;
    typedef value_type*                     pointer;
    typedef const value_type*               const_pointer;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef io::BasicBuffer                 buffer_type;

    // Function pointers
	typedef const U32 (self_type::*getStrideFunction)(const buffer_type& buffer, const U32 position) const;

public:
    InfoContainer();
    InfoContainer(const DataContainer& container);
    ~InfoContainer(void);

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
    inline reference back(void){ return(this->__containers[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__containers[this->n_entries - 1]); }

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__containers[0]); }
    inline iterator end()  { return iterator(&this->__containers[this->n_entries - 1]); }
    inline const_iterator begin()  const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->__containers[this->n_entries - 1]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->__containers[this->n_entries - 1]); }

private:
    template <class actual_primitive>
    void __setup(const DataContainer& container, getStrideFunction func){
		if(container.buffer_strides_uncompressed.size() == 0)
			return;

		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

		U32 current_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
			new( &this->__containers[i] ) value_type( container, current_offset, (this->*func)(container.buffer_strides_uncompressed, i) );
			current_offset += (this->*func)(container.buffer_strides_uncompressed, i) * sizeof(actual_primitive);
		}
		assert(current_offset == container.buffer_data_uncompressed.size());
	}

	template <class actual_primitive>
	void __setup(const DataContainer& container, const U32 stride_size){
		this->n_entries = container.buffer_data_uncompressed.size() / sizeof(actual_primitive);

		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

		U32 current_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
			new( &this->__containers[i] ) value_type( container, current_offset, stride_size );
			current_offset += stride_size * sizeof(actual_primitive);
		}
		assert(current_offset == container.buffer_data_uncompressed.size());
	}

	// Access function
	template <class stride_primitive> inline const U32 __getStride(const buffer_type& buffer, const U32 position) const{
		return(*reinterpret_cast<const stride_primitive* const>(&buffer.buffer[position*sizeof(stride_primitive)]));
	}

private:
    size_t  n_entries;
    pointer __containers;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
InfoContainer<return_type>::InfoContainer(const DataContainer& container) :
	n_entries(0),
	__containers(nullptr)
{
	if(container.buffer_data_uncompressed.size() == 0)
		return;

	if(container.header.controller.mixedStride){
		getStrideFunction func = nullptr;

		switch(container.header_stride.controller.type){
		case(tachyon::core::YON_TYPE_8B):  func = &self_type::__getStride<BYTE>; this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(BYTE); break;
		case(tachyon::core::YON_TYPE_16B): func = &self_type::__getStride<U16>;  this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U16);  break;
		case(tachyon::core::YON_TYPE_32B): func = &self_type::__getStride<U32>;  this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U32);  break;
		case(tachyon::core::YON_TYPE_64B): func = &self_type::__getStride<U64>;  this->n_entries = container.buffer_strides_uncompressed.size() / sizeof(U64);  break;
		default: std::cerr << "Disallowed stride" << std::endl; return;
		}

		if(container.header.isSigned()){
			switch(container.header.controller.type){
			case(tachyon::core::YON_TYPE_8B):  (this->__setup<SBYTE>(container, func)); break;
			case(tachyon::core::YON_TYPE_CHAR): (this->__setup<char>(container, func)); break;
			case(tachyon::core::YON_TYPE_16B): (this->__setup<S16>(container, func));   break;
			case(tachyon::core::YON_TYPE_32B): (this->__setup<S32>(container, func));   break;
			case(tachyon::core::YON_TYPE_64B): (this->__setup<S64>(container, func));   break;
			case(tachyon::core::YON_TYPE_FLOAT): (this->__setup<float>(container, func));   break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, func));   break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):  (this->__setup<BYTE>(container, func)); break;
			case(tachyon::core::YON_TYPE_16B): (this->__setup<U16>(container, func));  break;
			case(tachyon::core::YON_TYPE_32B): (this->__setup<U32>(container, func));  break;
			case(tachyon::core::YON_TYPE_64B): (this->__setup<U64>(container, func));  break;
			case(tachyon::core::YON_TYPE_FLOAT): (this->__setup<float>(container, func));   break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, func));   break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		}
	} else {
		if(container.header.controller.signedness){
			switch(container.header.controller.type){
			case(tachyon::core::YON_TYPE_8B):  (this->__setup<SBYTE>(container, container.header.stride)); break;
			case(tachyon::core::YON_TYPE_CHAR): (this->__setup<char>(container, container.header.stride)); break;
			case(tachyon::core::YON_TYPE_16B): (this->__setup<S16>(container, container.header.stride));   break;
			case(tachyon::core::YON_TYPE_32B): (this->__setup<S32>(container, container.header.stride));   break;
			case(tachyon::core::YON_TYPE_64B): (this->__setup<S64>(container, container.header.stride));   break;
			case(tachyon::core::YON_TYPE_FLOAT): (this->__setup<float>(container, container.header.stride));   break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, container.header.stride));   break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.controller.type){
			case(tachyon::core::YON_TYPE_8B):  (this->__setup<BYTE>(container, container.header.stride)); break;
			case(tachyon::core::YON_TYPE_16B): (this->__setup<U16>(container, container.header.stride));  break;
			case(tachyon::core::YON_TYPE_32B): (this->__setup<U32>(container, container.header.stride));  break;
			case(tachyon::core::YON_TYPE_64B): (this->__setup<U64>(container, container.header.stride));  break;
			case(tachyon::core::YON_TYPE_FLOAT): (this->__setup<float>(container, container.header.stride));   break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, container.header.stride));   break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		}
	}
}

template <class return_type>
InfoContainer<return_type>::~InfoContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__containers + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

}
}


#endif /* InfoContainer_H_ */
