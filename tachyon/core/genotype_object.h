#ifndef CORE_GENOTYPE_OBJECT_H_
#define CORE_GENOTYPE_OBJECT_H_

#include "meta_entry.h"

namespace tachyon{
namespace core{

#define YON_GT_RLE_ALLELE_A(PRIMITITVE, SHIFT, ADD)  (((PRIMITITVE) & ((1 << (SHIFT)) - 1) << (ADD)) >> (ADD));
#define YON_GT_RLE_ALLELE_B(PRIMITIVE, SHIFT, ADD)   (((PRIMITIVE) & ((1 << (SHIFT)) - 1) << ((ADD)+(SHIFT))) >> ((ADD)+(SHIFT)));
#define YON_GT_RLE_LENGTH(PRIMITIVE, SHIFT, ADD)     ((PRIMITIVE) >> (2*(SHIFT) + (ADD)))
#define YON_GT_DIPLOID_ALLELE_LOOKUP(A,B,shift,mask) (((A) & (mask)) << (shift)) | ((B) & (mask))
#define YON_GT_DIPLOID_BCF_A(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> ((SHIFT) + 1)) & ((1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_B(PRIMITIVE, SHIFT)       (((PRIMITIVE) >> 1) & ((1 << (SHIFT)) - 1))
#define YON_GT_DIPLOID_BCF_PHASE(PRIMITIVE)          ((PRIMITIVE) & 1)

const SBYTE YON_GT_RLE_CORRECTION[3] = {0, 0, 4};
const SBYTE YON_GT_RLE_RECODE[3] = {0,1,-2};

struct GTObjectAllele{
public:
	GTObjectAllele() : phase(0), allele(0){}
	~GTObjectAllele() = default;

public:
	BYTE  phase;
	SBYTE allele;
};

struct GTObject{
private:
    typedef GTObject             self_type;
    typedef GTObjectAllele       allele_type;
    typedef allele_type          value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::size_t          size_type;

public:
    GTObject(void) :
    	n_ploidy(0),
		n_objects(0),
		alleles(nullptr)
	{}

    GTObject(const self_type& other) :
    	n_ploidy(other.n_ploidy),
		n_objects(other.n_objects),
		alleles(new value_type[other.n_objects])
    {
		for(U32 i = 0; i < this->n_ploidy; ++i)
			this->alleles[i] = other.alleles[i];
    }

    GTObject(self_type&& other) noexcept :
    		n_ploidy(other.n_ploidy),
			n_objects(other.n_objects),
			alleles(other.alleles)
	{
		other.alleles = nullptr;
    }

    GTObject& operator=(const self_type& other){
		if(this->n_ploidy == other.n_ploidy){
			for(U32 i = 0; i < this->n_ploidy; ++i)
				this->alleles[i] = other.alleles[i];
		} else {
			delete [] this->alleles;
			this->alleles = new value_type[other.n_ploidy];

			for(U32 i = 0; i < other.n_ploidy; ++i)
				this->alleles[i] = other.alleles[i];
		}

		this->n_ploidy  = other.n_ploidy;
		this-> n_objects = other.n_objects;
		return(*this);
    }

    GTObject& operator=(self_type&& other) noexcept{
		if (this == &other)
			return *this;

		this->n_ploidy = other.n_ploidy;
		this->n_objects = other.n_objects;
		delete [] this->alleles;
		this->alleles = other.alleles;
		other.alleles = nullptr;
		return *this;
    }

    virtual ~GTObject(void){ delete [] this->alleles; }

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
    inline const size_t size(void) const{ return(this->n_objects); }

    // Element access
    inline pointer         data(void){ return(this->alleles); }
    inline const_pointer   data(void) const{ return(this->alleles); }
    inline reference       operator[](const U32& position){ return(this->alleles[position]); }
    inline const_reference operator[](const U32& position) const{ return(this->alleles[position]); }
    inline reference       at(const U32& position){ return(this->alleles[position]); }
    inline const_reference at(const U32& position) const{ return(this->alleles[position]); }

    // Iterator
    inline iterator       begin(){ return iterator(&this->alleles[0]); }
    inline iterator       end(){ return iterator(&this->alleles[this->n_objects]); }
    inline const_iterator begin() const{ return const_iterator(&this->alleles[0]); }
    inline const_iterator end() const{ return const_iterator(&this->alleles[this->n_objects]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->alleles[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->alleles[this->n_objects]); }

    friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
    	if(entry.n_ploidy){
    		if(entry.alleles[0].allele == -1){
    			stream.put('.');
    			return(stream);
    		}
    		if(entry.alleles[0].allele == -2) stream.put('.');
    		else stream << (int)entry.alleles[0].allele;

    		for(U32 i = 1; i < entry.n_ploidy; ++i){
    			if(entry.alleles[i].allele == -1) break;
    			stream << (entry.alleles[i].phase ? '|' : '/');
    			if(entry.alleles[i].allele == -2) stream.put('.');
    			else stream << (int)entry.alleles[i].allele;
    		}
    	}
    	return(stream);
    }

    friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& entry){
		if(entry.n_ploidy){
			if(entry.alleles[0].allele == -1){
				buffer += '.';
				return(buffer);
			}
			if(entry.alleles[0].allele == -2) buffer += '.';
			else buffer.AddReadble(entry.alleles[0].allele);

			for(U32 i = 1; i < entry.n_ploidy; ++i){
				if(entry.alleles[i].allele == -1) break;
				buffer += (entry.alleles[i].phase ? '|' : '/');
				if(entry.alleles[i].allele == -2) buffer += '.';
				else buffer.AddReadble(entry.alleles[i].allele);
			}
		}
		return(buffer);
	}

public:
    BYTE      n_ploidy; // ploidy of objects
    size_type n_objects;// number of objects
    pointer   alleles;  // alleleic data -> length = n_ploidy
};

struct GTObjectDiploidRLE : public GTObject{
private:
	typedef GTObjectDiploidRLE self_type;
	typedef core::MetaEntry    meta_type;

public:
    GTObjectDiploidRLE(void){}

    template <class T>
    void operator()(const T& gt_primitive, const meta_type& meta_entry)
    {
    	this->__interpret<T>(gt_primitive, meta_entry);
    }

    void operator()(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		this->__interpret(n_entries, alleleA, alleleB, phase);
	}

	void __interpret(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy  = 2;
			this->alleles    = new core::GTObjectAllele[2];
		}

		this->n_objects         = n_entries;
		this->alleles[0].allele  = alleleA;
		this->alleles[1].allele  = alleleB;
		this->alleles[0].phase = phase;
	}

private:
    template <class T>
    void __interpret(const T& gt_primitive, const meta_type& meta_entry)
    {
    	if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy  = 2;
			this->alleles   = new core::GTObjectAllele[2];
    	}
		const BYTE shift = meta_entry.isAnyGTMissing()   ? 2 : 1;
		const BYTE add   = meta_entry.isGTMixedPhasing() ? 1 : 0;

		if(add) this->alleles[0].phase = gt_primitive & 1;
		else    this->alleles[0].phase = meta_entry.getControllerPhase();

		this->n_objects        = YON_GT_RLE_LENGTH(gt_primitive, shift, add);
		this->alleles[0].allele = YON_GT_RLE_ALLELE_A(gt_primitive, shift, add);
		this->alleles[1].allele = YON_GT_RLE_ALLELE_B(gt_primitive, shift, add);
	}
};

struct GTObjectDiploidSimple : public GTObject{
private:
	typedef GTObjectDiploidSimple self_type;
	typedef core::MetaEntry       meta_type;

public:
	GTObjectDiploidSimple(void){}

    template <class T>
    void operator()(const T& gt_primitive, const meta_type& meta_entry)
    {
    	this->__interpret<T>(gt_primitive, meta_entry);
    }

    void operator()(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		this->__interpret(n_entries, alleleA, alleleB, phase);
	}

	void __interpret(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy  = 2;
			this->alleles   = new core::GTObjectAllele[2];
		}

		this->n_objects         = n_entries;
		this->alleles[0].allele  = alleleA;
		this->alleles[1].allele  = alleleB;
		this->alleles[0].phase = phase;
	}

private:
    template <class T>
    void __interpret(const T& gt_primitive, const meta_type& meta_entry)
    {
    	if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy     = 2;
			this->alleles      = new core::GTObjectAllele[2];
    	}
		const BYTE shift    = ceil(log2(meta_entry.getNumberAlleles() + 1 + meta_entry.isAnyGTMissing())); // Bits occupied per allele, 1 value for missing
		const BYTE add      = meta_entry.isGTMixedPhasing() ? 1 : 0;

		if(add) this->alleles[0].phase = gt_primitive & 1;
		else    this->alleles[0].phase = meta_entry.getControllerPhase();

		this->n_objects        = YON_GT_RLE_LENGTH(gt_primitive, shift, add);
		this->alleles[0].allele = YON_GT_RLE_ALLELE_A(gt_primitive, shift, add);
		this->alleles[1].allele = YON_GT_RLE_ALLELE_B(gt_primitive, shift, add);
	}
};

struct GTObjectDiploidBCF : public GTObject{
private:
	typedef GTObjectDiploidBCF self_type;
	typedef core::MetaEntry    meta_type;

public:
	GTObjectDiploidBCF(void){}

    template <class T>
    void operator()(const T& gt_primitive, const meta_type& meta_entry)
    {
    	this->__interpret<T>(gt_primitive, meta_entry);
    }

    void operator()(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		this->__interpret(n_entries, alleleA, alleleB, phase);
	}

	void __interpret(const U32& n_entries, const S32& alleleA, const S32& alleleB, const bool& phase)
	{
		if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy  = 2;
			this->alleles   = new core::GTObjectAllele[2];
		}

		this->n_objects         = n_entries;
		this->alleles[0].allele  = alleleA;
		this->alleles[1].allele  = alleleB;
		this->alleles[0].phase = phase;
	}

private:
    template <class T>
    void __interpret(const T& gt_primitive, const meta_type& meta_entry)
    {
    	if(this->alleles == nullptr || this->n_ploidy != 2){
			delete [] this->alleles;
			this->n_ploidy     = 2;
			this->alleles      = new core::GTObjectAllele[2];
    	}
		const BYTE shift    = (sizeof(T)*8 - 1) / 2;

		this->alleles[0].phase = YON_GT_DIPLOID_BCF_PHASE(gt_primitive);
		this->n_objects         = 1;
		this->alleles[0].allele  = YON_GT_DIPLOID_BCF_A(gt_primitive, shift);
		this->alleles[1].allele  = YON_GT_DIPLOID_BCF_B(gt_primitive, shift);
	}
};

}
}



#endif /* GTOBJECT_H_ */
