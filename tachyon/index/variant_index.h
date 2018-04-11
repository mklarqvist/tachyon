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
    VariantIndexBin() :
    	blockID(0),
		n_variants_(0),
		n_blocks_(0),
		n_capacity_(100),
		blocks_(new U32[this->capacity()])
	{

	}

    VariantIndexBin(const self_type& other);

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
    void Add(const U32& variant_block_number){
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

private:
    friend std::ostream& operator<<(std::ostream& stream, const self_type& bin){
		// Do not write if this bin is empty
    	if(bin.n_blocks_ == 0) return(stream);

    	stream.write(reinterpret_cast<const char*>(&bin.blockID), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&bin.n_variants_), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&bin.n_blocks_), sizeof(size_type));
		for(U32 i = 0; i < bin.size(); ++i)
			stream.write(reinterpret_cast<const char*>(&bin.blocks_[i]), sizeof(value_type));

		return(stream);
	}

public:
	U32       blockID;
	U32       n_variants_; // number of variants belonging to this bin
	size_type n_blocks_;
	size_type n_capacity_;
	pointer   blocks_;    // tachyon blocks belonging to this bin
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
    VariantIndexContig() :
    	l_contig_(0),
		l_contig_rounded_(0),
		n_bins_(0),
		n_capacity_(0),
		n_levels_(0),
		bins_cumsum_(nullptr),
		bins_(nullptr)
	{

	}

    VariantIndexContig(const U64 l_contig, const BYTE n_levels) :
    	l_contig_(l_contig),
		l_contig_rounded_(0),
		n_bins_(0),
		n_capacity_(0),
		n_levels_(n_levels),
		bins_cumsum_(nullptr),
		bins_(nullptr)
    {
    	this->l_contig_rounded_ = this->roundLengthClosestBase4_(this->l_contig_);
    	if(this->n_levels_ != 0){
    		this->calculateCumulativeSums_();
    		this->n_capacity_ = this->bins_cumsum_[this->n_levels_ - 1] + 64;
    		this->n_bins_ = this->bins_cumsum_[this->n_levels_ - 1];
    		this->bins_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
    		for(U32 i = 0; i < this->size(); ++i){
    			new( &this->bins_[i] ) value_type(  );
    			this->bins_[i].blockID = i;
    		}
    	}
    }

    VariantIndexContig(const self_type& other);

    ~VariantIndexContig(){
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->bins_ + i)->~VariantIndexBin();

		::operator delete[](static_cast<void*>(this->bins_));
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

	inline const S32 Add(const U64& fromPosition, const U64& toPosition, const U32& yon_block_id){
		for(S32 i = this->n_levels_ - 1; i != 0; --i){
			U32 binFrom = S64(fromPosition/(this->l_contig_rounded_ / pow(4,i)));
			U32 binTo = S64(toPosition/(this->l_contig_rounded_ / pow(4,i)));
			if(binFrom == binTo){
				//std::cerr << fromPosition << "->" << toPosition << ", adding to " << binFrom << " level " << i << " cum : " << this->bins_cumsum_[i-1]+binFrom << "/" << this->size() << std::endl;

				this->bins_[this->bins_cumsum_[i - 1]+binFrom].Add(yon_block_id);
				return(this->bins_cumsum_[i - 1]+binFrom);
			}
		}
		this->bins_[0].Add(yon_block_id);
		return(0);
	}

private:
    inline U64 roundLengthClosestBase4_(const U64& length) const{
		return( ( pow(4,this->n_levels_) - (length % (U64)pow(4,this->n_levels_)) ) + length );
    }

    void calculateCumulativeSums_(void){
    	if(this->n_levels_ == 0) return;

    	delete [] this->bins_cumsum_;
    	this->bins_cumsum_ = new U32[this->n_levels_];

    	U32 total = 0;
    	for(U32 i = 0; i < this->n_levels_; ++i){
    		total += pow(4,i);
    		this->bins_cumsum_[i] = total - 1; // remove 0 to start relative zero
    	}
    }

    friend std::ostream& operator<<(std::ostream& stream, const self_type& contig){
    	stream.write(reinterpret_cast<const char*>(&contig.l_contig_), sizeof(U64));
    	stream.write(reinterpret_cast<const char*>(&contig.l_contig_rounded_), sizeof(U64));
    	stream.write(reinterpret_cast<const char*>(&contig.n_bins_), sizeof(size_type));
    	stream.write(reinterpret_cast<const char*>(&contig.n_levels_), sizeof(BYTE));
    	for(U32 i = 0; i < contig.size(); ++i) stream << contig.bins_[i];
    	return(stream);
    }

private:
	U64       l_contig_;         // as described in header
	U64       l_contig_rounded_; // rounded up to next base-4
	size_type n_bins_;
	size_type n_capacity_;
	BYTE      n_levels_;    // 7 by default
	U32*      bins_cumsum_; // 1, 1+4, 1+4+16, 1+4+16+64, ...
	pointer   bins_;        // bin information
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
		contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{}

	VariantIndex(const self_type& other) :
		n_contigs_(other.n_contigs_),
		n_capacity_(other.n_capacity_),
		contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
	{
		for(U32 i = 0; i < this->size(); ++i)
			new( &this->contigs_[i] ) value_type( other.contigs_[i] );
	}

	~VariantIndex(){
		for(std::size_t i = 0; i < this->size(); ++i)
			(this->contigs_ + i)->~VariantIndexContig();

		::operator delete[](static_cast<void*>(this->contigs_));
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

	inline self_type& add(const U64& l_contig, const BYTE& n_levels){
		if(this->size() + 1 >= this->n_capacity_)
			this->resize();

		new( &this->contigs_[this->n_contigs_++] ) value_type( l_contig, n_levels );
		return(*this);
	}

	inline self_type& operator+=(const const_reference index_entry){
		if(this->size() + 1 >= this->n_capacity_)
			this->resize();

		new( &this->contigs_[this->n_contigs_++] ) value_type( index_entry );
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
	friend std::ostream& operator<<(std::ostream& stream, const self_type& index){
		stream.write(reinterpret_cast<const char*>(&index.n_contigs_), sizeof(size_type));
		for(U32 i = 0; i < index.size(); ++i) stream << index.contigs_[i];
		return(stream);
	}

private:
	size_type n_contigs_; // number of contigs
	size_type n_capacity_;
	pointer   contigs_;
};

}
}



#endif /* INDEX_VARIANT_INDEX_H_ */
