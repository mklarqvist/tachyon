#ifndef INDEX_VARIANT_INDEX_H_
#define INDEX_VARIANT_INDEX_H_

#include "variant_index_linear.h"

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
    	contigID_(0),
    	l_contig_(0),
		l_contig_rounded_(0),
		n_bins_(0),
		n_capacity_(0),
		n_levels_(0),
		n_sites_(0),
		bins_cumsum_(nullptr),
		bins_(nullptr)
	{

	}

    VariantIndexContig(const U32 contigID, const U64 l_contig, const BYTE n_levels) :
    	contigID_(contigID),
    	l_contig_(l_contig),
		l_contig_rounded_(0),
		n_bins_(0),
		n_capacity_(0),
		n_levels_(n_levels),
		n_sites_(0),
		bins_cumsum_(nullptr),
		bins_(nullptr)
    {
    	this->l_contig_rounded_ = this->roundLengthClosestBase4_(this->l_contig_);
    	if(this->n_levels_ != 0){
    		this->calculateCumulativeSums_();
    		this->n_capacity_ = this->bins_cumsum_[this->n_levels_] + 64;
    		this->n_bins_     = this->bins_cumsum_[this->n_levels_];
    		this->bins_       = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
    		for(U32 i = 0; i < this->size(); ++i){
    			new( &this->bins_[i] ) value_type(  );
    			this->bins_[i].binID_ = i;
    		}
    	}
    }

    VariantIndexContig(const self_type& other) :
    	contigID_(other.contigID_),
    	l_contig_(other.l_contig_),
		l_contig_rounded_(other.l_contig_rounded_),
		n_bins_(other.n_bins_),
		n_capacity_(other.n_capacity_),
		n_levels_(other.n_levels_),
		n_sites_(other.n_sites_),
		bins_cumsum_(nullptr),
		bins_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {
    	this->calculateCumulativeSums_();
    	for(U32 i = 0; i < this->size(); ++i)
    		new( &this->bins_[i] ) value_type( other.bins_[i] );
    }

    void operator=(const self_type& other){
    	// Clean previous
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->bins_ + i)->~VariantIndexBin();

		::operator delete[](static_cast<void*>(this->bins_));
		delete [] this->bins_cumsum_;

    	this->contigID_ = other.contigID_;
    	this->l_contig_ = other.l_contig_;
    	this->l_contig_rounded_ = other.l_contig_rounded_;
    	this->n_bins_ = other.n_bins_;
    	this->n_capacity_ = other.n_capacity_;
    	this->n_levels_ = other.n_levels_;
    	this->bins_cumsum_ = nullptr;
    	this->calculateCumulativeSums_();

    	this->bins_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
		for(U32 i = 0; i < this->size(); ++i)
			new( &this->bins_[i] ) value_type( other.bins_[i] );
    }


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
	inline const size_type& size_sites(void) const{ return(this->n_sites_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->bins_[0]); }
	inline iterator end(){ return iterator(&this->bins_[this->n_bins_]); }
	inline const_iterator begin() const{ return const_iterator(&this->bins_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->bins_[this->n_bins_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->bins_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->bins_[this->n_bins_]); }

	// Accessor
	inline U32& getContigID(void){ return(this->contigID_); }
	inline const U32& getContigID(void) const{ return(this->contigID_); }

	/**<
	 * Add a target interval tuple (from,to,block_ID)
	 * @param fromPosition From position of interval
	 * @param toPosition   To position of interval
	 * @param yon_block_id Tachyon block ID (generally a cumulative integer)
	 * @return
	 */
	inline const S32 Add(const U64& fromPosition, const U64& toPosition, const U32& yon_block_id){
		for(S32 i = this->n_levels_; i != 0; --i){
			U32 binFrom = S64(fromPosition/(this->l_contig_rounded_ / pow(4,i)));
			U32 binTo   = S64(toPosition/(this->l_contig_rounded_ / pow(4,i)));
			/**
			 * If both ends of the interval map into the same chunk we know the interval is
			 * completely contained: in this case we deposit the interval there
			 **/
			if(binFrom == binTo){
				//if(i != this->n_levels_) std::cerr << fromPosition << "->" << toPosition << ", adding to " << binFrom << " level " << i << " cum : " << this->bins_cumsum_[i-1]+binFrom << "/" << this->size() << std::endl;
				++this->n_sites_;
				this->bins_[this->bins_cumsum_[i - 1]+binFrom].Add(yon_block_id);
				return(this->bins_cumsum_[i - 1]+binFrom);
			}
		}
		this->bins_[0].Add(yon_block_id);
		++this->n_sites_;
		return(0);
	}

	/**<
	 * Computes the possible bins an interval might overlap
	 * @param from_position From position of interval
	 * @param to_position   To position of interval
	 * @return              Returns a vector of viable overlapping bins
	 */
	std::vector<value_type> possibleBins(const U64& from_position, const U64& to_position, const bool filter = true) const{
		std::vector<value_type> overlapping_chunks;
		//overlapping_chunks.push_back(this->at(0)); // level 0
		for(S32 i = this->n_levels_; i != 0; --i){
			U32 binFrom = S64(from_position/(this->l_contig_rounded_ / pow(4,i)));
			U32 binTo   = S64(to_position/(this->l_contig_rounded_ / pow(4,i)));

			//std::cerr << i << "/" << (int)this->n_levels_ << ": " << this->bins_cumsum_[i-1] << " + " << binFrom << "<>" << binTo << "/" << this->size() << std::endl;

			// Overlap from cumpos + (binFrom, binTo)
			// All these chunks could potentially hold intervals overlapping
			// the desired coordinates
			for(U32 j = binFrom; j <= binTo; ++j){
				if(filter == false)
					overlapping_chunks.push_back(this->at(this->bins_cumsum_[i - 1] + j));
				else {
					if(this->at(this->bins_cumsum_[i - 1] + j).size())
						overlapping_chunks.push_back(this->at(this->bins_cumsum_[i - 1] + j));
				}
			}
		}
		return(overlapping_chunks);
	}

private:
	/**<
	 * Round target integer up to the closest number divisible by 4
	 * @param length Input integer start value
	 * @return       Return a target integer divisible by 4
	 */
    inline U64 roundLengthClosestBase4_(const U64& length) const{
		return( ( pow(4,this->n_levels_) - (length % (U64)pow(4,this->n_levels_)) ) + length );
    }

    /**<
     * Pre-calculate the cumulative distribution of 4^(0:levels-1).
     * These values are used to find the array offset for levels > 0
     */
    void calculateCumulativeSums_(void){
    	if(this->n_levels_ == 0) return;

    	delete [] this->bins_cumsum_;
    	this->bins_cumsum_ = new U32[this->n_levels_ + 1]; // inclusive last

    	U32 total = 0;
    	for(U32 i = 0; i <= this->n_levels_; ++i){
    		total += pow(4,i);
    		this->bins_cumsum_[i] = total - 1; // remove 0 to start relative zero
    	}
    }

    friend std::ostream& operator<<(std::ostream& stream, const self_type& contig){
    	stream.write(reinterpret_cast<const char*>(&contig.contigID_),         sizeof(U32));
    	stream.write(reinterpret_cast<const char*>(&contig.l_contig_),         sizeof(U64));
    	stream.write(reinterpret_cast<const char*>(&contig.l_contig_rounded_), sizeof(U64));
    	stream.write(reinterpret_cast<const char*>(&contig.n_bins_),           sizeof(size_type));
    	stream.write(reinterpret_cast<const char*>(&contig.n_levels_),         sizeof(BYTE));
    	stream.write(reinterpret_cast<const char*>(&contig.n_sites_),          sizeof(size_type));

    	size_type n_items_written = 0;
    	for(U32 i = 0; i < contig.size(); ++i){
    		if(contig.bins_[i].size()) ++n_items_written;
    	}
    	stream.write(reinterpret_cast<const char*>(&n_items_written), sizeof(size_type));

    	for(U32 i = 0; i < contig.size(); ++i){
    		// If bins[i] contains data
    		if(contig.bins_[i].size())
    			stream << contig.bins_[i];
    	}
    	return(stream);
    }

    friend std::istream& operator>>(std::istream& stream, self_type& contig){
    	// Clear old data
    	if(contig.size()){
			for(std::size_t i = 0; i < contig.size(); ++i)
				(contig.bins_ + i)->~VariantIndexBin();

			::operator delete[](static_cast<void*>(contig.bins_));
    	}

		delete [] contig.bins_cumsum_;
		contig.bins_cumsum_ = nullptr;

		stream.read(reinterpret_cast<char*>(&contig.contigID_),         sizeof(U32));
		stream.read(reinterpret_cast<char*>(&contig.l_contig_),         sizeof(U64));
		stream.read(reinterpret_cast<char*>(&contig.l_contig_rounded_), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&contig.n_bins_),           sizeof(size_type));
		stream.read(reinterpret_cast<char*>(&contig.n_levels_),         sizeof(BYTE));
		stream.read(reinterpret_cast<char*>(&contig.n_sites_),          sizeof(size_type));
		contig.n_capacity_ = contig.n_bins_ + 64;
		size_type n_items_written = 0;
		stream.read(reinterpret_cast<char*>(&n_items_written), sizeof(size_type));

		// Allocate new
		contig.bins_ = static_cast<pointer>(::operator new[](contig.capacity()*sizeof(value_type)));
		for(U32 i = 0; i < contig.size(); ++i){
			new( &contig.bins_[i] ) value_type(  );
			contig.bins_[i].binID_ = i;
		}
		contig.calculateCumulativeSums_();

		// Load data accordingly
		for(U32 i = 0; i < n_items_written; ++i){
			value_type temp;
			stream >> temp;
			//std::cerr << "loading: " << temp.size() << " entries" << std::endl;
			contig.bins_[temp.binID_] = temp;
		}
		//std::cerr << std::endl;
		return(stream);
	}

private:
    U32       contigID_;
	U64       l_contig_;         // as described in header
	U64       l_contig_rounded_; // rounded up to next base-4
	size_type n_bins_;
	size_type n_capacity_;
	BYTE      n_levels_;    // 7 by default
	size_type n_sites_;
	U32*      bins_cumsum_; // 1, 1+4, 1+4+16, 1+4+16+64, ...
	pointer   bins_;        // bin information
};

class VariantIndex{
private:
	typedef VariantIndex       self_type;
    typedef std::size_t        size_type;
    typedef VariantIndexContig value_type;
    typedef VariantIndexLinear linear_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef IndexEntry         linear_entry_type;

public:
	VariantIndex() :
		n_contigs_(0),
		n_capacity_(1000),
		contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)))),
		linear_(static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type))))
	{

	}

	VariantIndex(const self_type& other) :
		n_contigs_(other.n_contigs_),
		n_capacity_(other.n_capacity_),
		contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)))),
		linear_(static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type))))
	{
		for(U32 i = 0; i < this->size(); ++i){
			new( &this->contigs_[i] ) value_type( other.contigs_[i] );
			new( &this->linear_[i] )  linear_type( other.linear_[i] );
		}
	}

	~VariantIndex(){
		for(std::size_t i = 0; i < this->size(); ++i){
			(this->contigs_ + i)->~VariantIndexContig();
			(this->linear_ + i)->~VariantIndexLinear();
		}

		::operator delete[](static_cast<void*>(this->contigs_));
		::operator delete[](static_cast<void*>(this->linear_));
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

	inline linear_type& linear_at(const size_type& contig_id){ return(this->linear_[contig_id]); }
	inline const linear_type& linear_at(const size_type& contig_id) const{ return(this->linear_[contig_id]); }

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

	/**<
	 * Add a contig with n_levels to the chain
	 * @param l_contig Length of contig
	 * @param n_levels Number of desired 4^N levels
	 * @return         Returns a reference of self
	 */
	inline self_type& add(const U32& contigID, const U64& l_contig, const BYTE& n_levels){
		if(this->size() + 1 >= this->n_capacity_)
			this->resize();

		new( &this->contigs_[this->n_contigs_] ) value_type( contigID, l_contig, n_levels );
		new( &this->linear_[this->n_contigs_] ) linear_type( contigID );
		++this->n_contigs_;
		return(*this);
	}

	/**<
	 * Add index entry to the linear index given a contig id
	 * @param contigID Target linear index at position contigID
	 * @param entry    Target index entry to push back onto the linear index vector
	 * @return         Returns a reference of self
	 */
	inline self_type& add(const U32& contigID, const linear_entry_type& entry){
		this->linear_[contigID] += entry;
		return(*this);
	}

	/**<
	 * Resizes the index to accept more contigs than currently allocated
	 * memory for. Resizes for the quad-tree index and the linear index
	 */
	void resize(void){
		pointer temp = this->contigs_;
		linear_type* temp_linear = this->linear_;

		this->n_capacity_ *= 2;
		this->contigs_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
		this->linear_  = static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type)));


		// Lift over values from old addresses
		for(U32 i = 0; i < this->size(); ++i){
			new( &this->contigs_[i] ) value_type( temp[i] );
			new( &this->linear_[i] )  linear_type( temp_linear[i] );
		}

		// Clear temp
		for(std::size_t i = 0; i < this->size(); ++i){
			(temp + i)->~VariantIndexContig();
			(temp_linear + i)->~VariantIndexLinear();
		}

		::operator delete[](static_cast<void*>(temp));
		::operator delete[](static_cast<void*>(temp_linear));
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& index){
		stream.write(reinterpret_cast<const char*>(&index.n_contigs_), sizeof(size_type));
		size_type n_items_written = 0;
		for(U32 i = 0; i < index.size(); ++i){
			if(index.contigs_[i].size_sites()) ++n_items_written;
		}
		stream.write(reinterpret_cast<const char*>(&n_items_written), sizeof(size_type));

		for(U32 i = 0; i < index.size(); ++i){
			// Write if contig[i] contains data
			if(index.contigs_[i].size_sites())
				stream << index.contigs_[i];
		}

		// Write linear index
		for(U32 i = 0; i < index.size(); ++i) stream << index.linear_[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& index){
		// Clear old data
		if(index.size()){
			for(std::size_t i = 0; i < index.size(); ++i){
				(index.contigs_ + i)->~VariantIndexContig();
				(index.linear_ + i)->~VariantIndexLinear();
			}

			::operator delete[](static_cast<void*>(index.contigs_));
			::operator delete[](static_cast<void*>(index.linear_));
		}

		stream.read(reinterpret_cast<char*>(&index.n_contigs_), sizeof(size_type));
		index.n_capacity_ = index.size() + 64;
		size_type n_items_written = 0;
		stream.read(reinterpret_cast<char*>(&n_items_written), sizeof(size_type));

		// Allocate new data
		index.contigs_ = static_cast<pointer>(::operator new[](index.capacity()*sizeof(value_type)));
		index.linear_ = static_cast<linear_type*>(::operator new[](index.capacity()*sizeof(linear_type)));
		for(U32 i = 0; i < index.size(); ++i) {
			new( &index.contigs_[i] ) value_type( );
			new( &index.linear_[i] ) linear_type( );
		}

		// Load data and update accordingly
		for(U32 i = 0; i < n_items_written; ++i){
			value_type temp;
			stream >> temp; // Read
			index.at(temp.getContigID()) = temp;
		}

		// Load linear index
		for(U32 i = 0; i < index.size(); ++i) stream >> index.linear_[i];

		return(stream);
	}

private:
	size_type    n_contigs_; // number of contigs
	size_type    n_capacity_;
	pointer      contigs_;
	linear_type* linear_;
};

}
}



#endif /* INDEX_VARIANT_INDEX_H_ */
