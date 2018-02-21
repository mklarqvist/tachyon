#ifndef CONTAINERS_INFO_CONTAINER_STRING_H_
#define CONTAINERS_INFO_CONTAINER_STRING_H_

#include "meta_container.h"
#include "info_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

/**<
 * InfoContainer for strings
 */
template <>
class InfoContainer<std::string>{
private:
    typedef InfoContainer        self_type;
    typedef std::string          value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef std::size_t          size_type;
    typedef io::BasicBuffer      buffer_type;
    typedef DataContainer        data_container_type;
    typedef MetaContainer        meta_container_type;
    typedef StrideContainer<U32> stride_container_type;

public:
    InfoContainer() :
    	n_entries(0),
		__containers(nullptr)
	{

	}

    InfoContainer(const data_container_type& container) :
    	n_entries(0),
		__containers(nullptr)
    {
    	if(container.header.hasMixedStride())
    		this->__setup(container);
    	else
    		this->__setup(container, container.header.stride);
    }

    InfoContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches) :
    	n_entries(0),
		__containers(nullptr)
    {
    	if(data_container.header.hasMixedStride())
			this->__setupBalanced(data_container, meta_container, pattern_matches);
		else
			this->__setupBalanced(data_container, meta_container, pattern_matches, data_container.header.stride);
    }

    ~InfoContainer(void){
    	for(std::size_t i = 0; i < this->n_entries; ++i)
    		((this->__containers + i)->~basic_string)();

    	::operator delete[](static_cast<void*>(this->__containers));
    }

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
    // For mixed strides
    void __setup(const data_container_type& container){
    	if(container.buffer_strides_uncompressed.size() == 0)
			return;

    	stride_container_type strides(container);
    	this->n_entries = strides.size();

		if(this->size() == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));


		U32 current_offset = 0;
		for(U32 i = 0; i < this->size(); ++i){
			new( &this->__containers[i] ) value_type(&container.buffer_data_uncompressed.data()[current_offset], strides[i]);
			current_offset += strides[i];
		}
		assert(current_offset == container.buffer_data_uncompressed.size());
    }

    void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches){
    	this->n_entries = meta_container.size();

		if(this->n_entries == 0)
			return;

		stride_container_type strides(data_container);
		this->__containers = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

		U32 current_offset = 0;
		U32 stride_offset = 0;
		for(U32 i = 0; i < this->size(); ++i){
			// Meta entry has no INFO
			if(meta_container[i].getInfoPatternID() == -1){
				new( &this->__containers[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].getInfoPatternID()]){
				new( &this->__containers[i] ) value_type(&data_container.buffer_data_uncompressed.data()[current_offset], strides[stride_offset]);
				current_offset += strides[stride_offset];
				++stride_offset;
			}
			// Otherwise place an empty
			else {
				new( &this->__containers[i] ) value_type( );
			}
		}
		assert(current_offset == data_container.buffer_data_uncompressed.size());
    }

    // For fixed strides
	void __setup(const data_container_type& container, const U32 stride_size){
		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

		U32 current_offset = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
			new( &this->__containers[i] ) value_type(&container.buffer_data_uncompressed.data()[current_offset], stride_size);
			current_offset += stride_size;
		}
		assert(current_offset == container.buffer_data_uncompressed.size());
    }

	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const U32 stride_size){
    	this->n_entries = meta_container.size();

		if(this->n_entries == 0)
			return;

		this->__containers = static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type)));

		if(data_container.header.isUniform() == false){
			U32 current_offset = 0;
			for(U32 i = 0; i < this->n_entries; ++i){
				// If there are no INFO fields
				if(meta_container[i].getInfoPatternID() == -1){
					new( &this->__containers[i] ) value_type( );
				} // If pattern matches
				else if(pattern_matches[meta_container[i].getInfoPatternID()]){
					new( &this->__containers[i] ) value_type(&data_container.buffer_data_uncompressed.data()[current_offset], stride_size);
					current_offset += stride_size;
				}
				// Otherwise place an empty
				else {
					new( &this->__containers[i] ) value_type( );
				}
			}
			assert(current_offset == data_container.buffer_data_uncompressed.size());
		}
		// Data is uniform
		else {
			for(U32 i = 0; i < this->n_entries; ++i){
				// If pattern matches
				if(pattern_matches[meta_container[i].getInfoPatternID()]){
					new( &this->__containers[i] ) value_type(data_container.buffer_data_uncompressed.data(), stride_size);
				}
				// Otherwise place an empty
				else {
					new( &this->__containers[i] ) value_type( );
				}
			}
		}
    }

private:
    size_t  n_entries;
    pointer __containers;
};

}
}



#endif /* CONTAINERS_INFO_CONTAINER_STRING_H_ */
