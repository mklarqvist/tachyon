#ifndef CORE_META_ALLELE_H_
#define CORE_META_ALLELE_H_

#include "io/basic_buffer.h"

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
	MetaAllele();
	MetaAllele(const self_type& other); // copy ctor
	MetaAllele(self_type&& other);
	MetaAllele(const char reference); // Ctor from packed byte
	MetaAllele(const std::string& reference);
	MetaAllele(const char* const in); // Ctor from buffer
	MetaAllele(const buffer_type& buffer, const U32 position); // Ctor directly from buffer object
	~MetaAllele(void);
	self_type& operator=(const self_type& other);

	void operator()(const char* const in);
	void operator()(const char* const in, const U32 length);

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
