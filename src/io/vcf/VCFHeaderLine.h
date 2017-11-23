#ifndef VCF_VCFHEADERLINE_H_
#define VCF_VCFHEADERLINE_H_

#include <algorithm> // for std::find

namespace Tachyon{
namespace VCF{

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
	VCFHeaderLine(const char* data, const U32 size) : size_(size), data_(data){}
	~VCFHeaderLine(){}

	inline const U32 size(void) const{ return this->pairs_.size(); }
	inline const key_value& operator[](const U32 p) const{ return this->pairs_[p]; }
	inline bool isValid(void) const{ return(this->size_ > 2 && (this->data_[0] == '#' && this->data_[1] == '#')); }
	inline bool isCONTIG(void) const{
		return(strncasecmp(&Constants::HEADER_CONTIG[0], &this->data_[0], Constants::HEADER_CONTIG.size()) == 0);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& pair){
		out << pair.data_ << '\n';
		for(U32 i = 0; i < pair.pairs_.size(); ++i)
			out << i << '/' << pair.size() << '\t' << pair[i] << '\n';

		return(out);
	}

	bool Parse(void){
		// Make sure this is a valid VCF header line
		// Rule: has to start with ##
		if(!this->isValid()){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Invalid VCF header line..." << std::endl;
			return false;
		}


		// Attempt to find an equal sign
		const char* match = std::find(this->data_, &this->data_[this->size_], '=');
		if(*match != '='){
			std::cerr << Helpers::timestamp("ERROR", "VCF") << "Corrupted VCF header entry: no equal match..." << std::endl;
			return false;
		}

		if(this->data_[match-this->data_+1] != '<'){
			std::cerr << "is text only: " << std::endl;
			std::cerr << std::string(this->data_, this->size_+1) << std::endl;
			return true;
		}

		// Find first equal sign
		// If next symbol is < then parse it
		std::string test(&this->data_[match-this->data_+2], this->size_ - (match-this->data_+2));
		std::vector<std::string> x = Helpers::split(test, ',');
		std::cerr << test << std::endl;

		for(U32 i = 0; i < x.size(); ++i){
			std::cerr << x[i] << std::endl;
			std::vector<std::string> y = Helpers::split(x[i], '=');

			for(U32 j = 0; j < y.size(); ++j)
				std::cerr << '\t' << y[j] << std::endl;

			if(y.size() != 2){
				std::cerr << "impossible: " << y.size() << std::endl;
				exit(1);
			}
			this->pairs_.push_back(key_value(y[0], y[1]));

		}
		return true;
	}

public:
	U32 size_;
	const char* data_;
	std::vector<key_value> pairs_;
};

}
}

#endif /* VCF_VCFHEADERLINE_H_ */
