#ifndef VCF_VCFHEADERLINE_H_
#define VCF_VCFHEADERLINE_H_

#include <iostream>
#include <cassert>
#include <algorithm> // for std::find
#include <vector>
#include <cstring>

#include "vcf_header_constants.h"
#include "support/helpers.h"

namespace tachyon{
namespace vcf{

enum TACHYON_VCF_HEADER_LINE_TYPE{
	YON_VCF_HEADER_UNKNOWN,
	YON_VCF_HEADER_FORMAT,
	YON_VCF_HEADER_FILTER,
	YON_VCF_HEADER_INFO,
	YON_VCF_HEADER_CONTIG
};

struct VCFHeaderLine{
private:
	typedef VCFHeaderLine self_type;

	// Internal helper struct
	struct VCFHeaderLineKeyValue{
		typedef VCFHeaderLineKeyValue self_type;

		VCFHeaderLineKeyValue(const std::string& key, const std::string& value) : KEY(key), VALUE(value){}
		VCFHeaderLineKeyValue(){}
		~VCFHeaderLineKeyValue(){}

		friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
			out << pair.KEY << '\t' << pair.VALUE;
			return(out);
		}

		std::string KEY;
		std::string VALUE;
	};
	typedef VCFHeaderLineKeyValue key_value;

public:
	VCFHeaderLine(const char* data, const U32 size);
	~VCFHeaderLine();

	inline U32 size(void) const{ return this->pairs.size(); }
	inline const key_value& operator[](const U32 p) const{ return this->pairs[p]; }
	inline bool isValid(void) const{ return(this->size_ > 2 && (this->data[0] == '#' && this->data[1] == '#')); }
	inline bool isCONTIG(void) const{
		return(strncasecmp(&constants::HEADER_CONTIG[0], &this->data[0], constants::HEADER_CONTIG.size()) == 0);
	}

	bool Parse(void);

private:
	bool nextKey(U32& startPos);

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		out << pair.data << '\n';
		for(U32 i = 0; i < pair.pairs.size(); ++i)
			out << i << '/' << pair.size() << '\t' << pair[i] << '\n';

		return(out);
	}

public:
	TACHYON_VCF_HEADER_LINE_TYPE type;
	bool isIndexable; // if this record should form part of the primary map
	U32 size_; // size in bytes
	const char* data; // pointer to data
	std::vector<key_value> pairs; // key-value pairs
};

}
}

#endif /* VCF_VCFHEADERLINE_H_ */
