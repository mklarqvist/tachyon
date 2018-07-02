#ifndef VCF_VCFHEADERLINE_H_
#define VCF_VCFHEADERLINE_H_

#include <iostream>
#include <cassert>
#include <algorithm> // for std::find
#include <vector>
#include <cstring>

#include "VCFHeaderConstants.h"
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
	VCFHeaderLine(const char* data, const U32 size) : type(YON_VCF_HEADER_UNKNOWN), isIndexable(false), size_(size), data(data){}
	~VCFHeaderLine(){}

	inline const U32 size(void) const{ return this->pairs.size(); }
	inline const key_value& operator[](const U32 p) const{ return this->pairs[p]; }
	inline bool isValid(void) const{ return(this->size_ > 2 && (this->data[0] == '#' && this->data[1] == '#')); }
	inline bool isCONTIG(void) const{
		return(strncasecmp(&constants::HEADER_CONTIG[0], &this->data[0], constants::HEADER_CONTIG.size()) == 0);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		out << pair.data << '\n';
		for(U32 i = 0; i < pair.pairs.size(); ++i)
			out << i << '/' << pair.size() << '\t' << pair[i] << '\n';

		return(out);
	}

	bool Parse(void){
		// Make sure this is a valid VCF header line
		// Rule: has to start with ##
		if(!this->isValid()){
			std::cerr << utility::timestamp("ERROR", "VCF") << "Invalid VCF header line..." << std::endl;
			return false;
		}

		// Attempt to find an equal sign
		const char* match = std::find(this->data, &this->data[this->size_], '=');
		if(*match != '='){
			std::cerr << utility::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match..." << std::endl;
			return false;
		}

		if(this->data[match - this->data + 1] != '<'){
			this->isIndexable = false;
			return true;
		}

		//std::string(this->data + 2, match - this->data - 2);
		if(strncasecmp(this->data + 2, "FORMAT", match - this->data - 2) == 0){
			this->type = YON_VCF_HEADER_FORMAT;
			this->isIndexable = true;
		} else if(strncasecmp(this->data + 2, "FILTER", match - this->data - 2) == 0){
			this->type = YON_VCF_HEADER_FILTER;
			this->isIndexable = true;
		} else if(strncasecmp(this->data + 2, "INFO", match - this->data - 2) == 0){
			this->type = YON_VCF_HEADER_INFO;
			this->isIndexable = true;
		}

		U32 matchPos = match - this->data + 1;
		if(this->data[matchPos] == '<'){
			if(this->data[this->size_] != '>'){
				std::cerr << utility::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: " << this->data[this->size_] << std::endl;
				return false;
			}

			++matchPos;

			// Sweep over and assert it is valid
			while(this->nextKey(matchPos)){
				// nothing in body
			}

		}

		return true;
	}

private:
	bool nextKey(U32& startPos){
		if(this->data[startPos] == '>')
			return false;

		const char* match = std::find(&this->data[startPos], &this->data[this->size_], '=');
		if(*match != '='){
			std::cerr << utility::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match in next key..." << std::endl;
			return false;
		}
		U32 matchPos = match - this->data;
		key_value entry;
		entry.KEY = std::string(&this->data[startPos], matchPos - startPos);

		startPos = matchPos + 1;

		char match_token = ',';
		BYTE adjust_value = 0;
		if(this->data[startPos] == '"'){
			match_token = '"';
			adjust_value = 1;
		}

		match = std::find(&this->data[startPos + adjust_value], &this->data[this->size_], match_token);
		if(*match == '>'){
			entry.VALUE = std::string(&this->data[startPos],this->size_ - startPos);
			startPos = matchPos + 1;
			this->pairs.push_back(entry);
			return false;
		} else if(*match != match_token){
			std::cerr << utility::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no comma match in next key..." << std::endl;
			return false;
		}

		matchPos = match - this->data;
		entry.VALUE = std::string(&this->data[startPos], matchPos - startPos + adjust_value);
		startPos = matchPos + 1;
		if(this->data[startPos] == ',') ++startPos;
		this->pairs.push_back(entry);
		return true;
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
