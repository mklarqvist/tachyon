#ifndef INDEX_VARIANT_INDEX_BIN_H_
#define INDEX_VARIANT_INDEX_BIN_H_

#include <cstring>
#include <fstream>

#include "support/type_definitions.h"

namespace tachyon{
namespace index{

struct VariantIndexBin{
private:
	typedef VariantIndexBin    self_type;
    typedef std::size_t        size_type;
    typedef U32                value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

public:
    VariantIndexBin() :
    	binID_(0),
		n_variants_(0),
		n_blocks_(0),
		n_capacity_(100),
		blocks_(new value_type[this->capacity()])
	{

	}

    VariantIndexBin(const self_type& other) :
    	binID_(other.binID_),
		n_variants_(other.n_variants_),
		n_blocks_(other.n_blocks_),
		n_capacity_(other.n_capacity_),
		blocks_(new value_type[this->capacity()])
    {
    	memcpy(this->blocks_, other.blocks_, sizeof(value_type)*other.n_blocks_);
    }


    VariantIndexBin& operator=(const self_type& other){
    	delete [] this->blocks_;
		this->blocks_     = new value_type[other.capacity()];
    	this->binID_      = other.binID_;
    	this->n_blocks_   = other.n_blocks_;
    	this->n_capacity_ = other.n_capacity_;
    	for(U32 i = 0; i < this->size(); ++i) this->blocks_[i] = other.blocks_[i];

    	return(*this);
    }

    inline bool operator<(const self_type& other) const{ return(this->binID_ < other.binID_); }

    ~VariantIndexBin(){ delete [] this->blocks_; }

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
	inline reference at(const size_type& position){ return(this->blocks_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->blocks_[position]); }
	inline reference operator[](const size_type& position){ return(this->blocks_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->blocks_[position]); }
	inline pointer data(void){ return(this->blocks_); }
	inline const_pointer data(void) const{ return(this->blocks_); }
	inline reference front(void){ return(this->blocks_[0]); }
	inline const_reference front(void) const{ return(this->blocks_[0]); }
	inline reference back(void){ return(this->blocks_[this->n_blocks_ - 1]); }
	inline const_reference back(void) const{ return(this->blocks_[this->n_blocks_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_blocks_ == 0); }
	inline const size_type& size(void) const{ return(this->n_blocks_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->blocks_[0]); }
	inline iterator end(){ return iterator(&this->blocks_[this->n_blocks_]); }
	inline const_iterator begin() const{ return const_iterator(&this->blocks_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->blocks_[this->n_blocks_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->blocks_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->blocks_[this->n_blocks_]); }

	// resize
	void resize(){
		pointer old = this->blocks_;
		this->n_capacity_ *= 2;
		this->blocks_ = new value_type[this->capacity()*2];
		for(U32 i = 0; i < this->size(); ++i) this->blocks_[i] = old[i];
		delete [] old;
	}

	/**<
	 * Update
	 * @param variant_block_number
	 */
    void add(const U32& variant_block_number){
		if(this->size() + 1 >= this->capacity())
			this->resize();

    	if(this->size()){ // Has data
    		if(this->back() != variant_block_number) // check parity between previous tachyon block and current one
    			this->blocks_[this->n_blocks_++] = variant_block_number;

    		++this->n_variants_;
    	} else { // Empty
    		this->blocks_[this->n_blocks_++] = variant_block_number;
    		++this->n_variants_;
    	}
    }

    std::ostream& print(std::ostream& stream){
    	stream << "ID: " << this->binID_ << ", variants: " << this->n_variants_ << ", associated blocks: " << this->n_blocks_;
    	if(this->size()){
    		stream << ", yon-blocks ids: " << this->blocks_[0];
    		for(U32 i = 1; i < this->size(); ++i)
    			stream << ',' << this->blocks_[i];
    	}

		return(stream);
	}

private:
    friend std::ostream& operator<<(std::ostream& stream, const self_type& bin){
		stream.write(reinterpret_cast<const char*>(&bin.binID_),     sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&bin.n_variants_), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&bin.n_blocks_),   sizeof(size_type));
		for(U32 i = 0; i < bin.size(); ++i)
			stream.write(reinterpret_cast<const char*>(&bin.blocks_[i]), sizeof(value_type));

		return(stream);
	}

    friend std::istream& operator>>(std::istream& stream, self_type& bin){
    	delete [] bin.blocks_;
 		stream.read(reinterpret_cast<char*>(&bin.binID_),     sizeof(U32));
		stream.read(reinterpret_cast<char*>(&bin.n_variants_), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&bin.n_blocks_),   sizeof(size_type));
		bin.n_capacity_ = bin.size() + 64;
		bin.blocks_ = new value_type[bin.capacity()];

		for(U32 i = 0; i < bin.size(); ++i)
			stream.read(reinterpret_cast<char*>(&bin.blocks_[i]), sizeof(value_type));

		return(stream);
	}

public:
	U32       binID_;
	U32       n_variants_; // number of variants belonging to this bin
	size_type n_blocks_;
	size_type n_capacity_;
	pointer   blocks_;    // tachyon blocks belonging to this bin
};

}
}



#endif /* INDEX_VARIANT_INDEX_BIN_H_ */
