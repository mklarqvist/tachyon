#ifndef INDEX_INDEXBLOCKENTRYBITVECTOR_H_
#define INDEX_INDEXBLOCKENTRYBITVECTOR_H_

namespace tachyon{
namespace containers{

// Size of entries in these records are
// inferred from the number of INFO/FORMAT/FILTER
// entries in all the records in a block
struct DataBlockBitvector{
	typedef DataBlockBitvector self_type;
    typedef std::size_t        size_type;
    typedef U32                value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;

public:
	DataBlockBitvector() : n_keys(0), local_keys(nullptr), bit_bytes(nullptr){}
	~DataBlockBitvector(){
		delete [] this->local_keys;
		delete [] this->bit_bytes;
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

	inline void update(const BYTE& value, const U32& pos){ this->bit_bytes[pos] = value; }

	inline void allocate(const U32& n_keys, const U32& n_bytes){
		delete [] this->local_keys;
		delete [] this->bit_bytes;
		this->local_keys = new U32[n_keys];
		this->bit_bytes = new BYTE[n_bytes];
		memset(this->bit_bytes, 0, n_bytes);
		this->n_keys = n_keys;
	}

	inline void allocate(const U32& n_bytes){
		delete [] this->bit_bytes;
		this->bit_bytes = new BYTE[n_bytes];
		memset(this->bit_bytes, 0, n_bytes);
	}

    // Element access
    inline reference key_at(const size_type& position){ return(this->local_keys[position]); }
    inline const_reference key_at(const size_type& position) const{ return(this->local_keys[position]); }
    inline pointer key_data(void){ return(this->local_keys); }
    inline const_pointer key_data(void) const{ return(this->local_keys); }
    inline reference key_front(void){ return(this->local_keys[0]); }
    inline const_reference key_front(void) const{ return(this->local_keys[0]); }
    inline reference key_back(void){ return(this->local_keys[this->n_keys - 1]); }
    inline const_reference key_back(void) const{ return(this->local_keys[this->n_keys - 1]); }

    // Bit access
    inline const bool operator[](const U32 position) const{ return((this->bit_bytes[position / 8] & (1 << (position % 8))) >> (position % 8)); }

    // Capacity
    inline const bool empty(void) const{ return(this->n_keys == 0); }
    inline const value_type& size(void) const{ return(this->n_keys); }

    // Iterator
    inline iterator begin(){ return iterator(&this->local_keys[0]); }
    inline iterator end()  { return iterator(&this->local_keys[this->n_keys - 1]); }
    inline const_iterator begin()  const{ return const_iterator(&this->local_keys[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->local_keys[this->n_keys]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->local_keys[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->local_keys[this->n_keys]); }

	// Utility
	inline const U32 getBaseSize(void) const{ return(sizeof(U32) + sizeof(U32)*this->n_keys); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.n_keys), sizeof(U32));
		for(U32 i = 0; i < entry.n_keys; ++i)
			stream.write(reinterpret_cast<const char*>(&entry.local_keys[i]), sizeof(U32));

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.n_keys), sizeof(U32));
		entry.local_keys = new U32[entry.n_keys];
		for(U32 i = 0; i < entry.n_keys; ++i)
			stream.read(reinterpret_cast<char*>(&entry.local_keys[i]), sizeof(U32));

		return(stream);
	}

public:
	value_type n_keys;
	pointer    local_keys;
	BYTE*      bit_bytes;
};

}
}

#endif /* INDEX_INDEXBLOCKENTRYBITVECTOR_H_ */
