#ifndef INDEX_VARIANT_INDEX_H_
#define INDEX_VARIANT_INDEX_H_

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

    void Add(const U32& variant_block_number);

public:
	U32     blockID;
	U32     n_entries_; // number of variants belonging to this bin
	U32     n_blocks_;
	U32     n_capacity_;
	pointer blocks_;
};

class VariantIndexContig{
private:
	typedef VariantIndexContig self_type;
    typedef std::size_t        size_type;
    typedef VariantIndexBin    value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

public:

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
	inline reference at(const size_type& position){ return(this->bins_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->bins_[position]); }
	inline reference operator[](const size_type& position){ return(this->bins_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->bins_[position]); }
	inline pointer data(void){ return(this->bins_); }
	inline const_pointer data(void) const{ return(this->bins_); }
	inline reference front(void){ return(this->bins_[0]); }
	inline const_reference front(void) const{ return(this->bins_[0]); }
	inline reference back(void){ return(this->bins_[this->n_bins_ - 1]); }
	inline const_reference back(void) const{ return(this->bins_[this->n_bins_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_bins_ == 0); }
	inline const size_type& size(void) const{ return(this->n_bins_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->bins_[0]); }
	inline iterator end(){ return iterator(&this->bins_[this->n_bins_]); }
	inline const_iterator begin() const{ return const_iterator(&this->bins_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->bins_[this->n_bins_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->bins_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->bins_[this->n_bins_]); }

    void Add(const U32& target_bin, const U32& variant_block_number);

private:
	U64     l_contig_; // as described in header
	U64     l_contig_rounded_; // rounded up to next base-4
	U32     n_bins_;
	U32     n_capacity_;
	BYTE    n_levels_;
	pointer bins_;
};

class VariantIndex{
private:
	typedef VariantIndex       self_type;
    typedef std::size_t        size_type;
    typedef VariantIndexContig value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

public:

	VariantIndex() :
		n_contigs_(0),
		n_capacity_(1000),
		contigs_(new value_type[this->n_capacity_])
	{}

	VariantIndex(const self_type& other) :
		n_contigs_(other.n_contigs_),
		n_capacity_(other.n_capacity_),
		contigs_(new value_type[this->n_capacity_])
	{
		for(U32 i = 0; i < this->size(); ++i) this->contigs_[i] = other.contigs_[i];
	}

	~VariantIndex(){
		delete [] this->contigs_;
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
	inline reference at(const size_type& position){ return(this->contigs_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->contigs_[position]); }
	inline reference operator[](const size_type& position){ return(this->contigs_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->contigs_[position]); }
	inline pointer data(void){ return(this->contigs_); }
	inline const_pointer data(void) const{ return(this->contigs_); }
	inline reference front(void){ return(this->contigs_[0]); }
	inline const_reference front(void) const{ return(this->contigs_[0]); }
	inline reference back(void){ return(this->contigs_[this->n_contigs_ - 1]); }
	inline const_reference back(void) const{ return(this->contigs_[this->n_contigs_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_contigs_ == 0); }
	inline const size_type& size(void) const{ return(this->n_contigs_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->contigs_[0]); }
	inline iterator end(){ return iterator(&this->contigs_[this->n_contigs_]); }
	inline const_iterator begin() const{ return const_iterator(&this->contigs_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->contigs_[this->n_contigs_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->contigs_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->contigs_[this->n_contigs_]); }

	inline self_type& operator+=(const const_reference index_entry){
		if(this->size() + 1 == this->n_capacity_)
			this->resize();


		this->contigs_[this->n_contigs_++] = index_entry;
		return(*this);
	}
	inline self_type& add(const const_reference index_entry){ return(*this += index_entry); }

	void resize(void){
		pointer temp = this->contigs_;

		this->n_capacity_ *= 2;
		this->contigs_ = new value_type[this->capacity()];

		// Lift over values from old addresses
		for(U32 i = 0; i < this->size(); ++i)
			this->contigs_[i] = temp[i];

		delete [] temp;
	}

private:
	size_type n_contigs_; // number of contigs
	size_type n_capacity_;
	pointer   contigs_;
};

}
}



#endif /* INDEX_VARIANT_INDEX_H_ */
