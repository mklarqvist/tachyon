#ifndef VCF_VCFHEADER_H_
#define VCF_VCFHEADER_H_

#include <algorithm>

#include "../../support/helpers.h"
#include "../reader.h"
#include "VCFHeaderConstants.h"
#include "VCFHeaderLine.h"
#include "../BasicBuffer.h"
#include "../../algorithm/OpenHashTable.h"
#include "../../core/base/header/HeaderContig.h"
#include "../../core/base/header/HeaderMapEntry.h"
#include "../../core/base/header/HeaderSample.h"

namespace Tachyon {
namespace VCF{

class VCFHeader {
	typedef VCFHeader self_type;
	typedef Hash::HashTable<std::string, S32> hash_table_type;
	typedef Hash::HashTable<S32, U32> hash_table_map_type;
	typedef Core::HeaderContig contig_type;
	typedef IO::BasicBuffer buffer_type;
	typedef VCFHeaderLine header_line_type;
	typedef Core::HeaderMapEntry map_entry_type;
	typedef Core::HeaderSample header_sample_type;

	enum VCF_ERROR_TYPE {VCF_PASS, VCF_ERROR_LINE1, VCF_ERROR_LINES, VCF_ERROR_SAMPLE, STREAM_BAD};

public:
	VCFHeader();
	~VCFHeader();

	void unsetBorrowedPointers(void){
		this->contigsHashTable = nullptr;
		this->sampleHashTable = nullptr;
	}

	inline bool good(void) const{ return(this->error_bit == VCF_PASS); }
	inline bool valid(void) const{ return(this->version > 0); }
	inline void setVersion(float version){ this->version = version; }
	inline const float& getVersion(void) const{ return(this->version); }
	inline U32 getContigs(void) const{ return this->contigs.size(); }
	//inline const contig_type& operator[](const U32 p) const{ return(this->contigs[p]); }
	inline contig_type& getContig(const U32 p){ return this->contigs[p]; }
	inline U32 getLines(void) const{ return this->lines.size(); }
	inline const U64& size(void) const{ return this->samples; }

	inline bool getContig(const std::string& contig, S32*& retValue) const{
		return(this->contigsHashTable->GetItem(&contig[0], &contig, retValue, contig.size()));
	}

	inline bool getSample(const std::string& sample, S32*& retValue) const{
		return(this->sampleHashTable->GetItem(&sample[0], &sample, retValue, sample.size()));
	}

	bool parse(reader& stream);
	bool parse(const char* const data, const U32& length);

	FORCE_INLINE const map_entry_type& operator[](const U32& p) const{ return(this->map[p]); }

private:
	// These functions are unsafe as they require contigHashTable to be
	// set prior to calling
	// no tests are made to check
	inline void addContig(const std::string& contig, U32 value){
		this->contigsHashTable->SetItem(&contig[0], &contig, value, contig.size());
	}

	inline void addSample(const std::string& sample){
		this->sampleNames.push_back(header_sample_type(sample));
		this->sampleHashTable->SetItem(&sample[0], &sample, this->sampleNames.size()-1, sample.size());
	}

	// Internal overload helpers
	bool checkLine(const char* data, const U32 length);
	bool buildContigTable(void);
	void buildSampleTable(U64 samples);
	bool parseFirstLine(reader& stream);
	bool parseHeaderLines(reader& stream);
	bool parseSampleLine(reader& stream);
	bool parseFirstLine(const char* const data, U32& offset);
	bool parseHeaderLines(const char* const data, U32& offset);
	bool parseSampleLine(const char* const data, U32& offset, const U32& length);

private:
	// Output only
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const U32 n_contigs = entry.contigs.size();
		const U32 n_samples = entry.sampleNames.size();
		const U32 n_meta_fields = entry.map.size();
		stream.write(reinterpret_cast<const char*>(&n_contigs),sizeof(U32));
		for(U32 i = 0; i < n_contigs; ++i)
			stream << entry.contigs[i];

		stream.write(reinterpret_cast<const char*>(&n_samples),sizeof(U32));
		for(U32 i = 0; i < n_samples; ++i)
			stream << entry.sampleNames[i];

		stream.write(reinterpret_cast<const char*>(&n_meta_fields),sizeof(U32));
		for(U32 i = 0; i < n_meta_fields; ++i)
			stream << entry.map[i];

		return(stream);
	}

public:
	VCF_ERROR_TYPE error_bit;               // parse error bit
	U64 samples;                            // number of samples
	float version;                          // VCF version
	std::string literal;                    // string copy of header data
	// Contigs are written to disk as an
	// object
	std::vector<contig_type> contigs;       // contigs
	// Sample names are written to disk as
	// U32 l_name; char[l_name]
	std::vector<header_sample_type> sampleNames;   // sample names
	std::vector<header_line_type> lines;    // header lines
	std::vector<std::string> literal_lines; // vcf line literals
	// These entries are read from disk as
	// U32 l_name; BYTE type; char[l_name]
	std::vector<map_entry_type> map;
	// Constructed during run-time
	hash_table_type* contigsHashTable;     // hash table for contig names
	hash_table_type* sampleHashTable;      // hash table for sample names
	hash_table_map_type* map_lookup;       // hash map from name to identifier
};

}
} /* namespace Tomahawk */

#endif /* VCFHEADER_H_ */
