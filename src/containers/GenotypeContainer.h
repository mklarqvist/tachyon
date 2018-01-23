#ifndef CONTAINERS_GENOTYPECONTAINER_H_
#define CONTAINERS_GENOTYPECONTAINER_H_

#include "Block.h"
#include "GenotypeContainerInterface.h"
#include "MetaContainer.h"

namespace Tachyon{
namespace Core{

class GenotypeContainer{
private:
    typedef GenotypeContainer             self_type;
    typedef GenotypeContainerInterface    value_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;
    typedef MetaContainer                 meta_container_type;
    typedef MetaEntry                     meta_type;
    typedef IO::BasicBuffer               buffer_type;

    // Function pointers
	typedef const U32 (self_type::*getNativeFuncDef)(const buffer_type& buffer, const U32 position) const;


public:
    GenotypeContainer(const Block& block) :
    	n_entries(0),
    	__meta_container(block),
    	__iterators(nullptr)
	{
		this->n_entries   = this->__meta_container.size();
		this->__iterators = static_cast<pointer>(::operator new[](this->n_entries * sizeof(value_type)));

		if(this->n_entries == 0){
			std::cerr << "no eentries" << std::endl;
			exit(1);
			return;
		}

		const char* const data_rle    = block.gt_rle_container.buffer_data_uncompressed.data;
		const char* const data_simple = block.gt_simple_container.buffer_data_uncompressed.data;
		//std::cerr << block.gt_support_data_container.getSizeUncompressed() << std::endl;
		//std::cerr << block.gt_support_data_container.header.controller.uniform << std::endl;
		//std::cerr << block.gt_support_data_container.header.controller.mixedStride << std::endl;
		//std::cerr << block.gt_support_data_container.buffer_strides_uncompressed.size() << std::endl;

		if(block.gt_support_data_container.buffer_data_uncompressed.size() == 0){
			std::cerr << "is 0" << std::endl;
			exit(1);
		}

		if(block.gt_support_data_container.buffer_strides_uncompressed.size() == 0){
			std::cerr << "stride is 0" << std::endl;
			exit(1);
		}

		// data (0: rle, 1: simple), strides (n_objects)
		getNativeFuncDef getTarget = nullptr;
		switch(block.gt_support_data_container.header.controller.type){
		case(YON_TYPE_8B):  getTarget = &self_type::getNative<BYTE>; break;
		case(YON_TYPE_16B): getTarget = &self_type::getNative<U16>; break;
		case(YON_TYPE_32B): getTarget = &self_type::getNative<U32>; break;
		case(YON_TYPE_64B): getTarget = &self_type::getNative<U64>; break;
		default: std::cerr << "illegal type" << std::endl; return;
		}

		getNativeFuncDef getObjects = nullptr;
		switch(block.gt_support_data_container.header_stride.controller.type){
		case(YON_TYPE_8B):  getObjects = &self_type::getNative<BYTE>; break;
		case(YON_TYPE_16B): getObjects = &self_type::getNative<U16>; break;
		case(YON_TYPE_32B): getObjects = &self_type::getNative<U32>; break;
		case(YON_TYPE_64B): getObjects = &self_type::getNative<U64>; break;
		default: std::cerr << "illegal type" << std::endl; return;
		}

		U32 current_offset_rle    = 0;
		U32 current_offset_simple = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			// new( &this->__iterators[i] ) value_type( &container.buffer_data_uncompressed.data[current_offset], getStride(i), this->__meta_container[i] );
			const U32 n_objects = (this->*getTarget)(block.gt_support_data_container.buffer_data_uncompressed, i);
			const U32 target    = (this->*getObjects)(block.gt_support_data_container.buffer_strides_uncompressed, i);
			//std::cerr << i << '/' << this->n_entries << '\t' << n_objects << '\t' << target << '\t' << (int)this->__meta_container[i].hot.getPrimitiveWidth() << std::endl;
			if(target == 1){
				new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				current_offset_rle += n_objects * this->__meta_container[i].hot.getPrimitiveWidth();
			} else if(target == 2){
				new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				current_offset_simple += n_objects * this->__meta_container[i].hot.getPrimitiveWidth();
			} else {
				std::cerr << "illegal" << std::endl;
				exit(1);
			}
		}

		//std::cerr << current_offset_rle << '/' << block.gt_rle_container.buffer_data_uncompressed.size() << std::endl;
		//std::cerr << current_offset_simple << '/' << block.gt_simple_container.buffer_data_uncompressed.size() << std::endl;
		assert(current_offset_rle == block.gt_rle_container.buffer_data_uncompressed.size());
		assert(current_offset_simple == block.gt_simple_container.buffer_data_uncompressed.size());
	}

    ~GenotypeContainer(){
    		for(std::size_t i = 0; i < this->n_entries; ++i)
    			(this->__iterators + i)->~GenotypeContainerInterface();

    		::operator delete[](static_cast<void*>(this->__iterators));
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

    // Capacity
	inline const bool       empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

	// Element access
	inline pointer         data(void){ return(this->__iterators); }
	inline const_pointer   data(void) const{ return(this->__iterators); }
	inline reference       operator[](const U32& position){ return(this->__iterators[position]); }
	inline const_reference operator[](const U32& position) const{ return(this->__iterators[position]); }
	inline reference       at(const U32& position){ return(this->__iterators[position]); }
	inline const_reference at(const U32& position) const{ return(this->__iterators[position]); }


private:
    template <class intrinsic_primitive> inline const U32 getNative(const buffer_type& buffer, const U32 position) const{
    		return(*reinterpret_cast<const intrinsic_primitive* const>(&buffer.data[position*sizeof(intrinsic_primitive)]));
    }

private:
    size_type           n_entries;
    meta_container_type __meta_container;
    pointer             __iterators;
};

}
}



#endif /* CONTAINERS_GENOTYPECONTAINER_H_ */