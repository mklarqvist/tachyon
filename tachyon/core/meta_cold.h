#ifndef CORE_BASE_META_COLD_H_
#define CORE_BASE_META_COLD_H_

#include <cassert>

#include "../io/bcf/BCFEntry.h"
#include "../containers/datacontainer.h"

namespace tachyon{
namespace core{

/** ColdMetaAllele:
 *  @brief Contains parts of the cold component of the hot-cold split of a variant site meta information
 *  This is a supportive structure. It keeps allele information
 *  as a typed string. This data structure is always cast
 *  directly from pre-loaded byte streams.
 */
struct ColdMetaAllele{
private:
	typedef ColdMetaAllele self_type;
	typedef io::BasicBuffer buffer_type;

public:
	ColdMetaAllele() :
		l_allele(0),
		allele(nullptr)
	{}

	~ColdMetaAllele(void){
		delete [] this->allele;
	}

	ColdMetaAllele(const self_type& other) :
		l_allele(other.l_allele),
		allele(new char[other.l_allele])
	{
		memcpy(this->allele, other.allele, other.l_allele);
	}

	ColdMetaAllele(ColdMetaAllele&& other) :
		l_allele(other.l_allele),
		allele(other.allele)
	{
		other.allele = nullptr;
	}

	ColdMetaAllele& operator=(const self_type& other){
		this->l_allele = other.l_allele;
		delete [] this->allele;
		this->allele = new char[other.l_allele];
		memcpy(this->allele, other.allele, other.l_allele);
		return *this;
	}

	void operator()(char* in){
		this->l_allele = *reinterpret_cast<U16*>(in);
		delete [] this->allele;
		this->allele   =  new char[this->l_allele];
		memcpy(this->allele, &in[sizeof(U16)], this->l_allele);
	}

	inline const U32 objectSize(void) const{ return(this->l_allele + sizeof(U16)); }
	inline const U16& size(void) const{ return(this->l_allele); }
	inline const std::string toString(void) const{ return(std::string(this->allele, this->l_allele)); }

private:
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		// Write out allele
		buffer += (U16)entry.l_allele;
		buffer.Add(entry.allele, entry.l_allele);

		return(buffer);
	}

public:
	U16 l_allele; /**< Byte length of allele data */
	char* allele; /**< Char array of allele */
};

// Do NOT reinterpret_cast this struct as an array
// as offsets needs to be interpreted
struct MetaCold{
private:
	typedef MetaCold self_type;
	typedef bcf::BCFEntry bcf_type;
	typedef containers::DataContainer stream_container;
	typedef ColdMetaAllele allele_type;

    typedef std::size_t       size_type;
    typedef allele_type       value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;


public:
	explicit MetaCold(void);
	MetaCold(char* in);
	~MetaCold(void);
	MetaCold(const self_type& other);
	MetaCold& operator=(const self_type& other);

	// Recycle or overload empty object
	void operator()(char* in);

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
	inline reference at(const size_type& position){ return(this->alleles[position]); }
	inline const_reference at(const size_type& position) const{ return(this->alleles[position]); }
	inline reference operator[](const size_type& position){ return(this->alleles[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->alleles[position]); }
	inline pointer data(void){ return(this->alleles); }
	inline const_pointer data(void) const{ return(this->alleles); }
	inline reference front(void){ return(this->alleles[0]); }
	inline const_reference front(void) const{ return(this->alleles[0]); }
	inline reference back(void){ return(this->alleles[this->n_allele - 1]); }
	inline const_reference back(void) const{ return(this->alleles[this->n_allele - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_allele == 0); }
	inline const U16& size(void) const{ return(this->n_allele); }

	//
	std::string getName(void) const{
		if(this->n_ID == 0)
			return(std::string("."));

		else return(std::string(this->ID, this->n_ID));
	}

	inline const char* getNamePointer(void) const{ return(this->ID); }
	inline const U16& getNameLength(void) const{ return(this->n_ID); }
	inline const U16& getNumberAlleles(void) const{ return(this->n_allele); }

	// Write out entry using BCF entry as template
	// and injects into buffer
	bool write(const bcf_type& entry, stream_container& buffer);

	std::vector<std::string> getAlleleStrings(void) const;

    // Iterator
    inline iterator begin(){ return iterator(&this->alleles[0]); }
    inline iterator end(){ return iterator(&this->alleles[this->n_allele - 1]); }
    inline const_iterator begin() const{ return const_iterator(&this->alleles[0]); }
    inline const_iterator end() const{ return const_iterator(&this->alleles[this->n_allele - 1]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->alleles[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->alleles[this->n_allele - 1]); }

public:
	/**< Byte length of this record. Is required in downstream iterators */
	U32 l_body;

	/**< Quality field */
	float QUAL;

	/**< Number of alleles */
	U16 n_allele;
	/**< Variant name (ID). Length is bytes
	 * and actual char* data
	 * Names are limited to 16 bits
	 */
	U16 n_ID;
	char* ID;

	// allele info
	// ALTs are limited to 16 bits each
	allele_type* alleles;
};

}
}

#endif /* CORE_BASE_META_COLD_H_ */
