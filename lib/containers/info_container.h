#ifndef CONTAINER_INFOCONTAINER_H_
#define CONTAINER_INFOCONTAINER_H_

#include "data_container.h"
#include "primitive_container.h"
#include "stride_container.h"
#include "meta_container.h"
#include "utility/support_vcf.h"
#include "info_container_interface.h"

namespace tachyon{
namespace containers{

template <class return_type>
class InfoContainer : public InfoContainerInterface{
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
    typedef DataContainer                   data_container_type;
    typedef MetaContainer                   meta_container_type;
    typedef StrideContainer<U32>            stride_container_type;

public:
    InfoContainer();
    InfoContainer(const bool is_flag);
    InfoContainer(const data_container_type& container);
    InfoContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches);
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

    // Iterator
    inline iterator begin(){ return iterator(&this->__containers[0]); }
    inline iterator end()  { return iterator(&this->__containers[this->n_entries]); }
    inline const_iterator begin()  const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->__containers[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->__containers[this->n_entries]); }

    // Type-specific
    inline std::ostream& to_vcf_string(std::ostream& stream, const U32 position) const{
    	//utility::to_vcf_string(stream, this->at(position).data(), this->at(position).size());
    	return(stream);
    }

    inline io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const U32 position) const{
    	utility::to_vcf_string(buffer, this->at(position).data(), this->at(position).size());
    	return(buffer);
    }

    inline io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const U32 position) const{
    	//utility::to_json_string(buffer, this->at(position));
		return(buffer);
    }

    inline bool emptyPosition(const U32& position) const{ return(this->at(position).empty()); }

private:
    // For mixed strides
    template <class actual_primitive>
    void __setup(const data_container_type& container);

    template <class actual_primitive>
	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches);

    // For fixed strides
	template <class actual_primitive>
	void __setup(const data_container_type& container, const U32 stride_size);

	template <class actual_primitive>
	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U32 stride_size);

private:
    pointer __containers;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
InfoContainer<return_type>::InfoContainer(void) :
	__containers(nullptr)
{

}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const bool is_flag) :
	__containers(static_cast<pointer>(::operator new[](1*sizeof(value_type))))
{
	this->n_entries  = 1;
	this->n_capacity = 1;
	// Set the primitive container value to 0. This
	// is required for the yon1_t structures to point
	// to something that is not simply a nullpointer.
	// It Has no other practical uses.
	new( &this->__containers[0] ) value_type( 0 );
}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const data_container_type& data_container,
                                          const meta_container_type& meta_container,
                                            const std::vector<bool>& pattern_matches) :
	__containers(nullptr)
{
	if(data_container.buffer_data_uncompressed.size() == 0){
		return;
	}

	if(data_container.header.data_header.HasMixedStride()){
		if(data_container.header.data_header.IsSigned()){
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setupBalanced<SBYTE>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_16B):    (this->__setupBalanced<S16>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_32B):    (this->__setupBalanced<S32>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_64B):    (this->__setupBalanced<S64>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setupBalanced<BYTE>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_16B):    (this->__setupBalanced<U16>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_32B):    (this->__setupBalanced<U32>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_64B):    (this->__setupBalanced<U64>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		}
	} else {
		if(data_container.header.data_header.IsSigned()){
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setupBalanced<SBYTE>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->__setupBalanced<S16>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_32B):    (this->__setupBalanced<S32>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_64B):    (this->__setupBalanced<S64>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setupBalanced<BYTE>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->__setupBalanced<U16>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_32B):    (this->__setupBalanced<U32>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_64B):    (this->__setupBalanced<U64>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		}
	}
}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const data_container_type& container) :
	__containers(nullptr)
{
	if(container.buffer_data_uncompressed.size() == 0)
		return;


	if(container.header.data_header.HasMixedStride()){
		if(container.header.data_header.IsSigned()){
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setup<SBYTE>(container));  break;
			case(YON_TYPE_16B):    (this->__setup<S16>(container));    break;
			case(YON_TYPE_32B):    (this->__setup<S32>(container));    break;
			case(YON_TYPE_64B):    (this->__setup<S64>(container));    break;
			case(YON_TYPE_FLOAT):  (this->__setup<float>(container));  break;
			case(YON_TYPE_DOUBLE): (this->__setup<double>(container)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setup<BYTE>(container));   break;
			case(YON_TYPE_16B):    (this->__setup<U16>(container));    break;
			case(YON_TYPE_32B):    (this->__setup<U32>(container));    break;
			case(YON_TYPE_64B):    (this->__setup<U64>(container));    break;
			case(YON_TYPE_FLOAT):  (this->__setup<float>(container));  break;
			case(YON_TYPE_DOUBLE): (this->__setup<double>(container)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		}
	} else {
		if(container.header.data_header.IsSigned()){
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setup<SBYTE>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->__setup<S16>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_32B):    (this->__setup<S32>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_64B):    (this->__setup<S64>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_FLOAT):  (this->__setup<float>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->__setup<double>(container, container.header.data_header.stride)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->__setup<BYTE>(container, container.header.data_header.stride));   break;
			case(YON_TYPE_16B):    (this->__setup<U16>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_32B):    (this->__setup<U32>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_64B):    (this->__setup<U64>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_FLOAT):  (this->__setup<float>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->__setup<double>(container, container.header.data_header.stride)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;

			}
		}
	}
}

template <class return_type>
InfoContainer<return_type>::~InfoContainer(){
	for(std::size_t i = 0; i < this->size(); ++i)
		((this->__containers + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::__setup(const data_container_type& container){
	if(container.buffer_strides_uncompressed.size() == 0)
		return;

	// Worst case number of entries
	this->n_capacity = container.buffer_data_uncompressed.size() / sizeof(actual_primitive);
	if(this->capacity() == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
	stride_container_type strides(container);

	U32 current_offset = 0;
	U32 i = 0;
	while(true){
		new( &this->__containers[i] ) value_type( container, current_offset, strides[i] );
		current_offset += strides[i] * sizeof(actual_primitive);
		++this->n_entries;
		++i;

		// Break condition
		if(current_offset == container.buffer_data_uncompressed.size())
			break;

		// Assertion of critical error
		assert(current_offset < container.buffer_data_uncompressed.size());
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::__setupBalanced(const data_container_type& data_container,
                                                 const meta_container_type& meta_container,
                                                   const std::vector<bool>& pattern_matches)
{
	this->n_entries = meta_container.size();
	if(this->size() == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));
	stride_container_type strides(data_container);

	U32 current_offset = 0;
	U32 stride_offset  = 0;

	for(U32 i = 0; i < this->size(); ++i){
		// There are no INFO fields
		if(meta_container[i].GetInfoPatternId() == -1){
			new( &this->__containers[i] ) value_type( );
		}
		// If pattern matches
		else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
			new( &this->__containers[i] ) value_type( data_container, current_offset, strides[stride_offset] );
			current_offset += strides[stride_offset] * sizeof(actual_primitive);
			++stride_offset;
		}
		// Otherwise place an empty
		else {
			new( &this->__containers[i] ) value_type( );
		}
	}

	assert(current_offset == data_container.buffer_data_uncompressed.size());
	assert(stride_offset == strides.size());
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::__setup(const data_container_type& container, const U32 stride_size){
	this->n_entries = container.buffer_data_uncompressed.size() / sizeof(actual_primitive);

	if(this->size() == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	U32 current_offset = 0;
	for(U32 i = 0; i < this->size(); ++i){
		//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
		new( &this->__containers[i] ) value_type( container, current_offset, stride_size );
		current_offset += stride_size * sizeof(actual_primitive);
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::__setupBalanced(const data_container_type& data_container,
                                                 const meta_container_type& meta_container,
                                                   const std::vector<bool>& pattern_matches,
                                                                 const U32  stride_size)
{
	this->n_entries = meta_container.size();
	if(this->size() == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	U32 current_offset = 0;
	// Case 1: if data is uniform
	if(data_container.header.data_header.IsUniform()){
		for(U32 i = 0; i < this->size(); ++i){
			// There are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->__containers[i] ) value_type( data_container, 0, stride_size );
			} else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		current_offset += stride_size * sizeof(actual_primitive);
	}
	// Case 2: if data is not uniform
	else {
		for(U32 i = 0; i < this->size(); ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, stride_size );
				current_offset += stride_size * sizeof(actual_primitive);
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
	}
	assert(current_offset == data_container.buffer_data_uncompressed.size());
}

}
}


#endif /* InfoContainer_H_ */
