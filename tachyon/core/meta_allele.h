#ifndef CORE_META_ALLELE_H_
#define CORE_META_ALLELE_H_

namespace tachyon{
namespace core{

/** MetaAllele:
 *  @brief Contains parts of the cold component of the hot-cold split of a variant site meta information
 *  This is a supportive structure. It keeps allele information
 *  as a typed string. This data structure is always cast
 *  directly from pre-loaded byte streams.
 */
struct MetaAllele{
private:
	typedef MetaAllele      self_type;
	typedef io::BasicBuffer buffer_type;

public:
	MetaAllele() :
		l_allele(0),
		allele(nullptr)
	{}

	// Ctor from packed byte
	MetaAllele(const char reference) :
		l_allele(1),
		allele(new char[this->l_allele])
	{
		this->allele[0] = reference;
	}

	MetaAllele(const std::string& reference) :
		l_allele(reference.size()),
		allele(new char[this->l_allele])
	{
		memcpy(this->allele, &reference[0], reference.size());
	}

	// Ctor from buffer
	MetaAllele(const char* const in) :
		l_allele(*reinterpret_cast<const U16* const>(in)),
		allele(new char[this->l_allele])
	{
		memcpy(this->allele, &in[sizeof(U16)], this->l_allele);
	}

	// Ctor directly from buffer object
	MetaAllele(const buffer_type& buffer, const U32 position) :
		l_allele(*reinterpret_cast<const U16* const>(&buffer[position])),
		allele(new char[this->l_allele])
	{
		memcpy(this->allele, &buffer[position + sizeof(U16)], this->l_allele);
	}

	~MetaAllele(void){
		delete [] this->allele;
	}

	MetaAllele(const self_type& other) :
		l_allele(other.l_allele),
		allele(new char[other.l_allele])
	{
		memcpy(this->allele, other.allele, other.l_allele);
	}

	MetaAllele(MetaAllele&& other) :
		l_allele(other.l_allele),
		allele(other.allele)
	{
		other.allele = nullptr;
	}

	self_type& operator=(const self_type& other){
		this->l_allele = other.l_allele;
		delete [] this->allele;
		this->allele = new char[other.l_allele];
		memcpy(this->allele, other.allele, other.l_allele);
		return *this;
	}

	void operator()(const char* const in){
		this->l_allele = *reinterpret_cast<const U16* const>(in);
		delete [] this->allele;
		this->allele   =  new char[this->l_allele];
		memcpy(this->allele, &in[sizeof(U16)], this->l_allele);
	}

	void operator()(const char* const in, const U32 length){
		this->l_allele = length;
		delete [] this->allele;
		this->allele   =  new char[this->l_allele];
		memcpy(this->allele, &in[0], this->l_allele);
	}

	inline const U16& size(void) const{ return(this->l_allele); }
	inline const U16& length(void) const{ return(this->l_allele); }
	inline const std::string toString(void) const{ return(std::string(this->allele, this->l_allele)); }

private:
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		// Write out allele
		buffer += (U16)entry.l_allele;
		buffer.Add(entry.allele, entry.l_allele);
		return(buffer);
	}

public:
	U16   l_allele;
	char* allele;
};

}
}




#endif /* CORE_META_ALLELE_H_ */
