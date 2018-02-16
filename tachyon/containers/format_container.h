#ifndef CONTAINERS_FORMAT_CONTAINER_H_
#define CONTAINERS_FORMAT_CONTAINER_H_

#include "datacontainer.h"
#include "primitive_group_container.h"
#include "meta_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

/**<
 * Primary class for FORMAT data
 */
template <class return_type>
class FormatContainer{
private:
    typedef FormatContainer                 self_type;
    typedef PrimitiveGroupContainer<return_type> value_type;
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
    FormatContainer();
    FormatContainer(const data_container_type& container, const U64 n_samples);
    FormatContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64 n_samples); // use when balancing
    ~FormatContainer(void);

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
    inline iterator end()  { return iterator(&this->__containers[this->n_entries]); }
    inline const_iterator begin()  const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->__containers[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->__containers[this->n_entries]); }

private:
    /**<
     *
     * @param container  Data container
     * @param n_samples  Number of samples
     * @param func       Function pointer to element accessor of stride data
     */
    template <class actual_primitive>
    void __setup(const data_container_type& container, const U64& n_samples);

    /**<
     *
     * @param data_container
     * @param meta_container
     * @param pattern_matches
     * @param n_samples
     * @param func
     */
    template <class actual_primitive>
	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples);

    /**<
     *
     * @param data_container
     * @param meta_container
     * @param pattern_matches
     * @param n_samples
     * @param stride_size
     */
    template <class actual_primitive>
    	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples, const U32 stride_size);

    /**<
     *
     * @param container   Data container
     * @param n_samples   Number of samples
     * @param stride_size Fixed stride size
     */
	template <class actual_primitive>
	void __setup(const data_container_type& container, const U64& n_samples, const U32 stride_size);

	// Access function
	template <class stride_primitive>
	inline const U32 __getStride(const buffer_type& buffer, const U32 position) const;

private:
    size_t  n_entries;
    pointer __containers;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
FormatContainer<return_type>::FormatContainer(const data_container_type& data_container,
                                              const meta_container_type& meta_container,
                                                const std::vector<bool>& pattern_matches,
                                                              const U64  n_samples) :
	n_entries(0),
	__containers(nullptr)
{
	if(data_container.buffer_data_uncompressed.size() == 0)
		return;

	if(data_container.header.hasMixedStride()){
		if(data_container.header.isSigned()){
			switch(data_container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setupBalanced<SBYTE>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setupBalanced<S16>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setupBalanced<S32>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setupBalanced<S64>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, n_samples));  break;
			default: std::cerr << "Disallowed type: " << (int)data_container.header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setupBalanced<BYTE>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setupBalanced<U16>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setupBalanced<U32>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setupBalanced<U64>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, n_samples));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, n_samples));  break;
			default: std::cerr << "Disallowed type: " << (int)data_container.header.controller.type << std::endl; return;
			}
		}
	} else {
		if(data_container.header.isSigned()){
			switch(data_container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setupBalanced<SBYTE>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setupBalanced<S16>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setupBalanced<S32>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setupBalanced<S64>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			default: std::cerr << "Disallowed type: " << (int)data_container.header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setupBalanced<BYTE>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setupBalanced<U16>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setupBalanced<U32>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setupBalanced<U64>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setupBalanced<float>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setupBalanced<double>(data_container, meta_container, pattern_matches, n_samples, data_container.header.stride));  break;
			default: std::cerr << "Disallowed type: " << (int)data_container.header.controller.type << std::endl; return;
			}
		}
	}
}

template <class return_type>
FormatContainer<return_type>::FormatContainer(const data_container_type& container, const U64 n_samples) :
	n_entries(0),
	__containers(nullptr)
{
	if(container.buffer_data_uncompressed.size() == 0)
		return;

	if(container.header.controller.mixedStride){
		if(container.header.isSigned()){
			switch(container.header.controller.type){
			case(tachyon::core::YON_TYPE_8B):     (this->__setup<SBYTE>(container, n_samples));  break;
			case(tachyon::core::YON_TYPE_CHAR):   (this->__setup<char>(container, n_samples));   break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setup<S16>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setup<S32>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setup<S64>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setup<float>(container, n_samples));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, n_samples)); break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setup<BYTE>(container, n_samples));   break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setup<U16>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setup<U32>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setup<U64>(container, n_samples));    break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setup<float>(container, n_samples));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, n_samples)); break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		}
	} else {
		if(container.header.isSigned()){
			switch(container.header.controller.type){
			case(tachyon::core::YON_TYPE_8B):     (this->__setup<SBYTE>(container, n_samples, container.header.getStride()));  break;
			case(tachyon::core::YON_TYPE_CHAR):   (this->__setup<char>(container, n_samples, container.header.getStride()));   break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setup<S16>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setup<S32>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setup<S64>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setup<float>(container, n_samples, container.header.getStride()));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, n_samples, container.header.getStride())); break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.getPrimitiveType()){
			case(tachyon::core::YON_TYPE_8B):     (this->__setup<BYTE>(container, n_samples, container.header.getStride()));   break;
			case(tachyon::core::YON_TYPE_16B):    (this->__setup<U16>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_32B):    (this->__setup<U32>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_64B):    (this->__setup<U64>(container, n_samples, container.header.getStride()));    break;
			case(tachyon::core::YON_TYPE_FLOAT):  (this->__setup<float>(container, n_samples, container.header.getStride()));  break;
			case(tachyon::core::YON_TYPE_DOUBLE): (this->__setup<double>(container, n_samples, container.header.getStride())); break;
			default: std::cerr << "Disallowed type: " << (int)container.header.controller.type << std::endl; return;
			}
		}
	}
}

template <class return_type>
FormatContainer<return_type>::~FormatContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__containers + i)->~PrimitiveGroupContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

template <class return_type>
template <class actual_primitive>
void FormatContainer<return_type>::__setup(const data_container_type& container, const U64& n_samples){
	if(container.buffer_strides_uncompressed.size() == 0)
		return;

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
	stride_container_type strides(container);

	U32 current_offset = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << (this->*func)(container.buffer_strides_uncompressed, i) << std::endl;
		new( &this->__containers[i] ) value_type( container, current_offset, n_samples, strides[i] );
		current_offset += strides[i] * sizeof(actual_primitive) * n_samples;
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void FormatContainer<return_type>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples){
		this->n_entries = meta_container.size();
		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));
		stride_container_type strides(data_container);

		U32 current_offset = 0;
		U32 strides_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].getFormatPatternID()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, strides[strides_offset] );
				current_offset += strides[strides_offset] * sizeof(actual_primitive);
				++strides_offset;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.buffer_data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void FormatContainer<return_type>::__setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U64& n_samples, const U32 stride_size){
		this->n_entries = meta_container.size();
	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

	U32 current_offset = 0;
	// Case 1: if data is uniform
	if(data_container.header.isUniform()){
		for(U32 i = 0; i < this->n_entries; ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].getFormatPatternID()]){
				new( &this->__containers[i] ) value_type( data_container, 0, n_samples, stride_size );
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}

		current_offset += stride_size * sizeof(actual_primitive);
	}
	// Case 2: if data is not uniform
	else {
		for(U32 i = 0; i < this->n_entries; ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].getFormatPatternID()]){
				new( &this->__containers[i] ) value_type( data_container, current_offset, n_samples, stride_size );
				current_offset += stride_size * sizeof(actual_primitive) * n_samples;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
	}
	assert(current_offset == data_container.buffer_data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void FormatContainer<return_type>::__setup(const data_container_type& container, const U64& n_samples, const U32 stride_size){
	this->n_entries = container.buffer_data_uncompressed.size() / sizeof(actual_primitive) / n_samples / stride_size;

	if(this->n_entries == 0)
		return;

	this->__containers = static_cast<pointer>(::operator new[](this->n_entries * sizeof(value_type)));

	U32 current_offset = 0;
	// Case 1: data is uniform -> give all samples the same value
	if(container.header.isUniform()){
		for(U32 i = 0; i < this->n_entries; ++i)
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );

	}
	// Case 2: data is not uniform -> interpret data
	else {
		for(U32 i = 0; i < this->n_entries; ++i){
			//std::cerr << current_offset << '/' << container.buffer_data_uncompressed.size() << '\t' << "fixed: " << stride_size << std::endl;
			new( &this->__containers[i] ) value_type( container, current_offset, n_samples, stride_size );
			current_offset += stride_size * sizeof(actual_primitive) * n_samples;
		}
	}
	assert(current_offset == container.buffer_data_uncompressed.size());
}

// Access function
template <class return_type>
template <class stride_primitive> inline const U32 FormatContainer<return_type>::__getStride(const buffer_type& buffer, const U32 position) const{
	return(*reinterpret_cast<const stride_primitive* const>(&buffer.buffer[position * sizeof(stride_primitive)]));
}

}
}


#endif /* CONTAINERS_FORMAT_CONTAINER_H_ */
