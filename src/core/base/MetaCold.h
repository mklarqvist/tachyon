#ifndef CORE_BASE_METACOLD_H_
#define CORE_BASE_METACOLD_H_

#include <cassert>

#include "../../io/bcf/BCFEntry.h"
#include "../StreamContainer.h"

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
	ColdMetaAllele() : l_allele(0), allele(nullptr){}
	ColdMetaAllele(char* in) :
		l_allele(*reinterpret_cast<U16*>(&in[0])),
		allele(&in[sizeof(U16)])
	{}
	~ColdMetaAllele(void){} // do nothing

	void operator()(char* in){
		this->l_allele = *reinterpret_cast<U16*>(in);
		this->allele = reinterpret_cast<char*>(&in[sizeof(U16)]);
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
	typedef Core::StreamContainer stream_container;
	typedef ColdMetaAllele allele_type;

public:
	explicit MetaCold(void);
	~MetaCold(void);

	void operator()(char* in){
		this->QUAL = *reinterpret_cast<float*>(in);
		this->n_allele = *reinterpret_cast<U16*>(&in[sizeof(float)]);
		this->n_ID = *reinterpret_cast<U16*>(&in[sizeof(float)+sizeof(U16)]);
		this->ID = &in[sizeof(float)+2*sizeof(U16)];
		this->alleles = new allele_type[this->n_allele];
		U32 cumpos = sizeof(float)+2*sizeof(U16)+this->n_ID;
		for(U32 i = 0; i < this->n_allele; ++i){
			this->alleles[i](&in[cumpos]);
			cumpos += this->alleles[i].objectSize();
		}
	}

	// Parse everything in this entry
	bool parse(void);

	// Parse the name only
	bool parseID(void);

	// Parse the alleles only
	bool parseAlleles(void);

	// Write out entry using BCF entry as template
	// and injects into buffer
	bool write(const bcf_type& entry, stream_container& buffer);

public:
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
