#ifndef VCF_VCFHEADER_H_
#define VCF_VCFHEADER_H_

#include <algorithm>
#include <unordered_map>

#include "algorithm/OpenHashTable.h"
#include "core/header/header_contig.h"
#include "core/header/header_map_entry.h"
#include "core/header/header_sample.h"
#include "io/basic_reader.h"
#include "vcf_header_constants.h"
#include "vcf_header_line.h"
#include "support/helpers.h"
#include "io/basic_buffer.h"

namespace tachyon {
namespace vcf{

class VCFHeader {
	typedef VCFHeader                         self_type;
	typedef hash::HashTable<std::string, S32> hash_table_type;
	typedef core::HeaderContig                contig_type;
	typedef io::BasicBuffer                   buffer_type;
	typedef VCFHeaderLine                     header_line_type;
	typedef core::HeaderMapEntry              map_entry_type;
	typedef core::HeaderSample                header_sample_type;
	typedef io::BasicReader                   reader_type;

	enum VCF_ERROR_TYPE {VCF_PASS, VCF_ERROR_LINE1, VCF_ERROR_LINES, VCF_ERROR_SAMPLE, STREAM_BAD};

public:
	VCFHeader();
	~VCFHeader();

	inline bool good(void) const{ return(this->error_bit == VCF_PASS); }
	inline bool valid(void) const{ return(this->vcf_version > 0); }
	inline void setVersion(float version){ this->vcf_version = version; }
	inline const float& getVersion(void) const{ return(this->vcf_version); }
	inline U32 getContigs(void) const{ return this->contigs.size(); }
	//inline const contig_type& operator[](const U32 p) const{ return(this->contigs[p]); }
	inline contig_type& getContig(const U32 p){ return this->contigs[p]; }
	inline U32 getLines(void) const{ return this->lines.size(); }
	inline const U64& size(void) const{ return this->n_samples; }

	inline bool getContig(const std::string& contig, S32& retValue) const{
		std::unordered_map<std::string,S32>::const_iterator ret = this->htable_contigs.find(contig);
		if(ret != this->htable_contigs.end()){
			retValue = ret->second;
			return(true);
		}
		return(false);
	}

	inline bool getSample(const std::string& sample, S32& retValue) const{
		std::unordered_map<std::string,S32>::const_iterator ret = this->htable_samples.find(sample);
		if(ret != this->htable_samples.end()){
			retValue = ret->second;
			return(true);
		}
		return(false);
	}

	bool parse(reader_type& stream);
	bool parse(const char* const data, const U32& length);

private:
	inline void addContig(const std::string& contig, const U32& value){
		if(this->contigs.size() == 0){
			this->htable_contigs[contig] = value;
			return;
		}

		if(this->htable_contigs.find(contig) == this->htable_contigs.end()){
			this->htable_contigs[contig] = value;
		}
	}

	inline void addSample(const std::string& sample){
		if(this->n_samples == 0){
			this->sample_names.push_back(header_sample_type(sample));
			this->htable_samples[sample] = this->n_samples;
			++this->n_samples;
			return;
		}

		if(this->htable_samples.find(sample) == this->htable_samples.end()){
			this->sample_names.push_back(header_sample_type(sample));
			this->htable_samples[sample] = this->n_samples;
			++this->n_samples;
		}
	}

	// Internal overload helpers
	bool checkLine(const char* data, const U32 length);
	bool parseFirstLine(reader_type& stream);
	bool parseHeaderLines(reader_type& stream);
	bool parseSampleLine(reader_type& stream);
	bool parseFirstLine(const char* const data, U32& offset);
	bool parseHeaderLines(const char* const data, U32& offset);
	bool parseSampleLine(const char* const data, U32& offset, const U32& length);

public:
	VCF_ERROR_TYPE error_bit;               // parse error bit
	U64 n_samples;                          // number of samples
	float vcf_version;                      // VCF version
	std::string literal;                    // string copy of header data
	// Contigs are written to disk as an
	// object
	std::vector<contig_type> contigs;       // contigs
	// Sample names are written to disk as
	// U32 l_name; char[l_name]
	std::vector<header_sample_type> sample_names;   // sample names
	std::vector<header_line_type> lines;    // header lines
	std::vector<std::string> literal_lines; // vcf line literals
	// These entries are read from disk as
	// U32 l_name; BYTE type; char[l_name]
	//std::vector<map_entry_type> map;
	std::vector<map_entry_type> info_map;
	std::vector<map_entry_type> format_map;
	std::vector<map_entry_type> filter_map;
	// Constructed during run-time
	S32* info_remap; // map from IDX to local id for O(1) lookup
	S32* format_remap;
	S32* filter_remap;

	std::unordered_map<std::string, S32> htable_contigs;
	std::unordered_map<std::string, S32> htable_samples;
	//hash_table_type* contigsHashTable;     // hash table for contig names
	//hash_table_type* sampleHashTable;      // hash table for sample names
};

}
} /* namespace Tomahawk */

#endif /* VCFHEADER_H_ */
