#ifndef CONTAINERS_GENOTYPECONTAINERINTERFACE_H_
#define CONTAINERS_GENOTYPECONTAINERINTERFACE_H_

#include "Container.h"

namespace Tachyon{
namespace Core{

class GenotypeContainerInterface{
private:
    typedef GenotypeContainerInterface  self_type;
    typedef std::size_t                 size_type;

public:
    GenotypeContainerInterface(void) : n_entries(0), __data(nullptr){}
    GenotypeContainerInterface(const char* const data, const size_type& n_entries, const U32& n_bytes) : n_entries(n_entries), __data(new char[n_bytes]){ memcpy(this->__data, data, n_bytes); }
    virtual ~GenotypeContainerInterface(){ delete [] this->__data; }

    // GT-specific functionality
    //virtual void getGTSummary(void) const =0;
    //virtual void getGTSummaryGroups(void) const =0;
    //virtual void std::vector<float> getAlleleFrequency(void) =0;
    //virtual void std::vector<bool> getSamplesMissingness(void) =0;
    //virtual void std::vector<U32> getSamplesPloidy(void) =0;
    //virtual void std::vector<sample_summary> getSamplesSummary(void) =0;
    //virtual void std::vector<upp_triagonal> compareSamplesPairwise(void) =0;

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

protected:
    size_type n_entries;
    char* __data;    // GT primitive
};

template <class T>
class GenotypeContainerDiploidRLE : public GenotypeContainerInterface{
private:
	typedef GenotypeContainerInterface    parent_type;
    typedef GenotypeContainerDiploidRLE   self_type;
    typedef T                             value_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;
    typedef std::ptrdiff_t                difference_type;
    typedef std::size_t                   size_type;
    typedef MetaEntry                     meta_type;

public:
    GenotypeContainerDiploidRLE() : __local(nullptr){}
    GenotypeContainerDiploidRLE(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
    		parent_type(data, n_entries, n_entries*sizeof(value_type)),
		__local(reinterpret_cast<const T* const>(this->__data))
	{

	}
    ~GenotypeContainerDiploidRLE(){ }

    void operator()(const char* const data, const U32 n_entries, const meta_type& meta_entry){
    		this->n_entries = n_entries;
    		delete [] this->__data;

    		const T* const re = reinterpret_cast<const T* const>(data);
		for(U32 i = 0; i < n_entries; ++i)
			this->__data[i] = re[i];
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
    inline reference at(const size_type& position){ return(this->__local[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__local[position]); }
    inline reference operator[](const size_type& position){ return(this->__local[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__local[position]); }
    inline pointer data(void){ return(this->__local); }
    inline const_pointer data(void) const{ return(this->__local); }
    inline reference front(void){ return(this->__local[0]); }
    inline const_reference front(void) const{ return(this->__local[0]); }
    inline reference back(void){ return(this->__local[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__local[this->n_entries - 1]); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__local[0]); }
    inline iterator end(){ return iterator(&this->__local[this->n_entries - 1]); }
    inline const_iterator begin() const{ return const_iterator(&this->__local[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__local[this->n_entries - 1]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__local[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__local[this->n_entries - 1]); }

private:
    const_pointer __local;
};

template <class T>
class GenotypeContainerDiploidSimple : public GenotypeContainerInterface{
private:
	typedef GenotypeContainerInterface     parent_type;
    typedef GenotypeContainerDiploidSimple self_type;
    typedef T                              value_type;
    typedef value_type&                    reference;
    typedef const value_type&              const_reference;
    typedef value_type*                    pointer;
    typedef const value_type*              const_pointer;
    typedef std::ptrdiff_t                 difference_type;
    typedef std::size_t                    size_type;
    typedef MetaEntry                      meta_type;

public:
    GenotypeContainerDiploidSimple() : __local(nullptr){}
    GenotypeContainerDiploidSimple(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
		parent_type(data, n_entries, n_entries*sizeof(value_type)),
		__local(reinterpret_cast<const T* const>(this->__data))
    	{

    	}
    ~GenotypeContainerDiploidSimple(){  }

    void operator()(const char* const data, const U32 n_entries, const meta_type& meta_entry){
    		this->n_entries = n_entries;
    		delete [] this->__data;

    		const T* const re = reinterpret_cast<const T* const>(data);
		for(U32 i = 0; i < n_entries; ++i)
			this->__data[i] = re[i];
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
    inline reference at(const size_type& position){ return(this->__local[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__local[position]); }
    inline reference operator[](const size_type& position){ return(this->__local[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__local[position]); }
    inline pointer data(void){ return(this->__local); }
    inline const_pointer data(void) const{ return(this->__local); }
    inline reference front(void){ return(this->__local[0]); }
    inline const_reference front(void) const{ return(this->__local[0]); }
    inline reference back(void){ return(this->__local[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__local[this->n_entries - 1]); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__local[0]); }
    inline iterator end(){ return iterator(&this->__local[this->n_entries - 1]); }
    inline const_iterator begin() const{ return const_iterator(&this->__local[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__local[this->n_entries - 1]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__local[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__local[this->n_entries - 1]); }

private:
    const_pointer __local;
};

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
				//this->__iterators[i] = new GenotypeContainerDiploidRLE<BYTE>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				current_offset_rle += n_objects * this->__meta_container[i].hot.getPrimitiveWidth();
			} else if(target == 2){
				new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				//this->__iterators[i] = new GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
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
    		for(std::size_t i = 0; i < this->n_entries; ++i){
    			(this->__iterators + i)->~GenotypeContainerInterface();
    		}
    		::operator delete[](static_cast<void*>(this->__iterators));
    }

    // Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

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

#endif /* CONTAINERS_GENOTYPECONTAINERINTERFACE_H_ */
