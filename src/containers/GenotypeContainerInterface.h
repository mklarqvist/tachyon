#ifndef CONTAINERS_GENOTYPECONTAINERINTERFACE_H_
#define CONTAINERS_GENOTYPECONTAINERINTERFACE_H_

#include "Container.h"
#include "../core/GTObject.h"

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
    virtual U32 getSum(void) const =0;

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
    typedef MetaHotController             hot_controller_type;

public:
    GenotypeContainerDiploidRLE() : __local(nullptr){}
    GenotypeContainerDiploidRLE(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
    		parent_type(data, n_entries, n_entries*sizeof(value_type)),
		__local(reinterpret_cast<const T* const>(this->__data)),
		__controller(meta_entry.hot.controller)
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

    // GT-specific
    U32 getSum(void) const{
    		//std::cerr << (void*)this->__data << '\t' << (void*)this->__local << std::endl;
    		U32 count = 0;
    		//std::cerr << "before accessing controller" << std::endl;
    		const BYTE shift    = this->__controller.gt_anyMissing    ? 2 : 1;
    		const BYTE add      = this->__controller.gt_mixed_phasing ? 1 : 0;

    		for(U32 i = 0; i < this->n_entries; ++i){
    			//std::cerr << i << "/" << this->n_entries << std::endl;
			count += this->at(i) >> (2*shift + add);
    		}
    		return(count);
    }

private:
    const_pointer __local;
    hot_controller_type __controller;
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
    typedef MetaHotController              hot_controller_type;

public:
    GenotypeContainerDiploidSimple() : n_alleles(0), __local(nullptr){}
    GenotypeContainerDiploidSimple(const char* const data, const U32 n_entries, const meta_type& meta_entry) :
    		parent_type(data, n_entries, n_entries*sizeof(value_type)),
		n_alleles(meta_entry.cold.n_allele),
		__local(reinterpret_cast<const T* const>(this->__data)),
		__controller(meta_entry.hot.controller)
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

    // GT-specific
	U32 getSum(void) const{
		/*
		U32 count = 0;

		const BYTE shift    = ceil(log2(this->n_alleles + this->__controller.gt_anyMissing)); // Bits occupied per allele, 1 value for missing
		const BYTE add      = this->__controller.gt_mixed_phasing ? 1 : 0;

		for(U32 i = 0; i < this->n_entries; ++i){
			count += this->at(i) >> (2*shift + add);
		}
		return(count);
		*/
		return(0);
	}

private:
	BYTE n_alleles;
    const_pointer __local;
    hot_controller_type __controller;
};


}
}

#endif /* CONTAINERS_GENOTYPECONTAINERINTERFACE_H_ */
