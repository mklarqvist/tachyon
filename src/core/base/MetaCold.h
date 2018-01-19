#ifndef CORE_BASE_METACOLD_H_
#define CORE_BASE_METACOLD_H_

#include <cassert>

#include "../../io/bcf/BCFEntry.h"
#include "../Container.h"

namespace Tachyon{
namespace Core{

/** ColdMetaAllele:
 *  @brief Contains parts of the cold component of the hot-cold split of a variant site meta information
 *  This is a supportive structure. It keeps allele information
 *  as a typed string. This data structure is always cast
 *  directly from pre-loaded byte streams.
 */
struct ColdMetaAllele{
private:
	typedef ColdMetaAllele self_type;
	typedef IO::BasicBuffer buffer_type;

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
	typedef BCF::BCFEntry bcf_type;
	typedef Core::Container stream_container;
	typedef ColdMetaAllele allele_type;

public:
	explicit MetaCold(void);
	MetaCold(char* in);
	~MetaCold(void);
	MetaCold(const self_type& other);
	MetaCold& operator=(const self_type& other);

	// Recycle or overload empty object
	void operator()(char* in);

	// Write out entry using BCF entry as template
	// and injects into buffer
	bool write(const bcf_type& entry, stream_container& buffer);

	std::vector<std::string> getAlleleStrings(void) const;

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

#endif /* CORE_BASE_METACOLD_H_ */
