#ifndef INDEX_VARIANT_INDEX_BIN_H_
#define INDEX_VARIANT_INDEX_BIN_H_

#include <cstring>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "generic_iterator.h"

namespace tachyon{
namespace index{

struct VariantIndexBin{
public:
	typedef VariantIndexBin    self_type;
    typedef std::size_t        size_type;
    typedef uint32_t           value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

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


    VariantIndexBin& operator=(const self_type& other) {
    	delete [] this->blocks_;
		this->blocks_     = new value_type[other.capacity()];
    	this->binID_      = other.binID_;
    	this->n_blocks_   = other.n_blocks_;
    	this->n_capacity_ = other.n_capacity_;
    	for (uint32_t i = 0; i < this->size(); ++i) this->blocks_[i] = other.blocks_[i];

    	return(*this);
    }

    ~VariantIndexBin() { delete [] this->blocks_; }

    self_type& operator+=(const self_type& other) {
    	this->n_variants_ += other.n_variants_;
    	if (this->n_blocks_ + other.n_blocks_ > this->capacity())
    		this->resize(this->n_blocks_ + other.n_blocks_ + 64);

    	for (int i = 0; i < other.n_blocks_; ++i) {
    		//std::cerr << "adding block " << other.blocks_[i] << " to offset " << this->n_ << std::endl;
    		this->blocks_[this->n_blocks_ + i] = other.blocks_[i];
    	}

    	this->n_blocks_ += other.n_blocks_;

    	// Sort and dedupe.
    	if (this->size() > 1) {
			std::sort(&this->blocks_[0], &this->blocks_[this->n_blocks_]);

			value_type* temp = new value_type[this->size()];
			temp[0] = this->blocks_[0];
			int new_size = 1;
			for (int i = 1; i < this->size(); ++i) {
				if (temp[new_size - 1] != this->blocks_[i]) {
					temp[new_size++] = this->blocks_[i];
				}
			}
			memcpy(temp, this->blocks_, new_size*sizeof(value_type));
			delete [] temp;
			this->n_blocks_ = new_size;
    	}

    	return(*this);
    }

    inline bool operator<(const self_type& other) const { return(this->binID_ < other.binID_); }

	// Element access
	inline reference at(const size_type& position) { return(this->blocks_[position]); }
	inline const_reference at(const size_type& position) const { return(this->blocks_[position]); }
	inline reference operator[](const size_type& position) { return(this->blocks_[position]); }
	inline const_reference operator[](const size_type& position) const { return(this->blocks_[position]); }
	inline pointer data(void) { return(this->blocks_); }
	inline const_pointer data(void) const { return(this->blocks_); }
	inline reference front(void) { return(this->blocks_[0]); }
	inline const_reference front(void) const { return(this->blocks_[0]); }
	inline reference back(void) { return(this->blocks_[this->n_blocks_ - 1]); }
	inline const_reference back(void) const { return(this->blocks_[this->n_blocks_ - 1]); }

	// Capacity
	inline bool empty(void) const { return(this->n_blocks_ == 0); }
	inline const size_type& size(void) const { return(this->n_blocks_); }
	inline const size_type& capacity(void) const { return(this->n_capacity_); }

	// Iterator
	inline iterator begin() { return iterator(&this->blocks_[0]); }
	inline iterator end() { return iterator(&this->blocks_[this->n_blocks_]); }
	inline const_iterator begin() const { return const_iterator(&this->blocks_[0]); }
	inline const_iterator end() const { return const_iterator(&this->blocks_[this->n_blocks_]); }
	inline const_iterator cbegin() const { return const_iterator(&this->blocks_[0]); }
	inline const_iterator cend() const { return const_iterator(&this->blocks_[this->n_blocks_]); }

	// resize
	void resize() {
		pointer old = this->blocks_;
		this->n_capacity_ *= 2;
		this->blocks_ = new value_type[this->capacity()*2];
		for (uint32_t i = 0; i < this->size(); ++i) this->blocks_[i] = old[i];
		delete [] old;
	}

	void resize(const uint32_t new_size) {
		pointer old = this->blocks_;
		this->n_capacity_ = new_size;
		this->blocks_ = new value_type[this->capacity()*2];
		for (uint32_t i = 0; i < this->size(); ++i) this->blocks_[i] = old[i];
		delete [] old;
	}

	/**<
	 * Update
	 * @param variant_block_number
	 */
    void Add(const uint32_t& variant_block_number) {
		if (this->size() + 1 >= this->capacity())
			this->resize();

    	if (this->size()) { // Has data
    		if (this->back() != variant_block_number) // check parity between previous tachyon block and current one
    			this->blocks_[this->n_blocks_++] = variant_block_number;

    		++this->n_variants_;
    	} else { // Empty
    		this->blocks_[this->n_blocks_++] = variant_block_number;
    		++this->n_variants_;
    	}
    }

    std::ostream& Print(std::ostream& stream) {
    	stream << "ID: " << this->binID_ << ", variants: " << this->n_variants_ << ", associated blocks: " << this->n_blocks_;
    	if (this->size()) {
    		stream << ", yon-blocks ids: " << this->blocks_[0];
    		for (uint32_t i = 1; i < this->size(); ++i)
    			stream << ',' << this->blocks_[i];
    	}

		return(stream);
	}

private:
    friend std::ostream& operator<<(std::ostream& stream, const self_type& bin) {
		stream.write(reinterpret_cast<const char*>(&bin.binID_),     sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&bin.n_variants_), sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&bin.n_blocks_),   sizeof(size_type));
		for (uint32_t i = 0; i < bin.size(); ++i)
			stream.write(reinterpret_cast<const char*>(&bin.blocks_[i]), sizeof(value_type));

		return(stream);
	}

    friend std::istream& operator>>(std::istream& stream, self_type& bin) {
    	delete [] bin.blocks_;
 		stream.read(reinterpret_cast<char*>(&bin.binID_),     sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&bin.n_variants_), sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&bin.n_blocks_),   sizeof(size_type));
		bin.n_capacity_ = bin.size() + 64;
		bin.blocks_ = new value_type[bin.capacity()];

		for (uint32_t i = 0; i < bin.size(); ++i)
			stream.read(reinterpret_cast<char*>(&bin.blocks_[i]), sizeof(value_type));

		return(stream);
	}

public:
	uint32_t  binID_;
	uint32_t  n_variants_; // number of variants belonging to this bin
	size_type n_blocks_;
	size_type n_capacity_;
	pointer   blocks_;    // tachyon blocks belonging to this bin
};

}
}



#endif /* INDEX_VARIANT_INDEX_BIN_H_ */
