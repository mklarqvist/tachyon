/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "support/magic_constants.h"
#include "io/basic_buffer.h"

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

#define YON_BYTE_MISSING  INT8_MIN
#define YON_BYTE_EOV      (INT8_MIN+1)
#define YON_SHORT_MISSING INT16_MIN
#define YON_SHORT_EOV     (INT16_MIN+1)
#define YON_INT_MISSING   INT32_MIN
#define YON_INT_EOV       (INT32_MIN+1)
#define YON_FLOAT_NAN     0x7FC00000
#define YON_FLOAT_MISSING 0x7F800001
#define YON_FLOAT_EOV     0x7F800002

namespace tachyon {

struct VcfContig {
public:
	VcfContig() : idx(0), n_bases(0){}
	~VcfContig() = default;

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##contig=<ID=GL000241.1,assembly=b37,length=42152>
		std::string ret = "##contig=<ID=" + this->name;
		if(extra.size()){
			ret += "," + this->extra[0].first + "=" + this->extra[0].second;
			for(uint32_t i = 1; i < this->extra.size(); ++i){
				ret += "," + this->extra[i].first + "=" + this->extra[i].second;
			}
		}
		if(this->description.size()) ret += ",Description=" + this->description;
		ret += ",length=" + std::to_string(this->n_bases);
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The name of the contig. Canonically this is the first
	// non-whitespace-containing string after the > marker in a FASTA file.
	// For example, the line:
	//      >chr1 more info here
	// has a name of "chr1" and a description of "more info here"
	std::string name;

	// Ideally this record is filled in as described above, but not all FASTA
	// readers capture the description information after the name. Since a
	// description is not required by the FASTA spec, we cannot distinguish cases
	// where a description was not present and where a parser ignored it.
	std::string description;

	// The length of this contig in basepairs.
	int64_t n_bases;

	// Additional information used when reading and writing VCF headers. An
	// example map of key-value extra fields would transform an input line
	// containing 'assembly=B36,taxonomy=x,species="Homo sapiens"' to a map with
	// "assembly" -> "B36", "taxonomy" -> "x", "species" -> "Homo sapiens". We
	// never use this information internally, other than reading it in so we can
	// write the contig out again.
	std::vector< std::pair<std::string, std::string> > extra;
};

// Temp declare
struct VcfInfo{
public:
	VcfInfo() : idx(0){}
	~VcfInfo() = default;

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
		std::string ret = "##INFO=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(this->source.size()) ret += ",Source=" + this->source;
		if(this->source.size()) ret += ",Version=" + this->version;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
		std::string ret = "##INFO=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(this->source.size()) ret += ",Source=" + this->source;
		if(this->source.size()) ret += ",Version=" + this->version;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the INFO field. Examples include "MQ0" or "END".
	std::string id;

	// Required. The number of values included with the info field. This should be
	// the string representation of the number, e.g. "1" for a single entry, "2"
	// for a pair of entries, etc. Special cases arise when the number of entries
	// depend on attributes of the Variant or are unknown in advance, and include:
	// "A": The field has one value per alternate allele.
	// "R": The field has one value per allele (including the reference).
	// "G": The field has one value for each possible genotype.
	// ".": The number of values varies, is unknown, or is unbounded.
	std::string number;

	// Required. The type of the INFO field. Valid values are "Integer", "Float",
	// "Flag", "Character", and "String".
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;

	// Optional. The annotation source used to generate the field.
	std::string source;

	// Optional. The version of the annotation source used to generate the field.
	std::string version;
};

struct VcfFormat{
public:
	VcfFormat() : idx(0){}
	~VcfFormat() = default;

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FORMAT=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FORMAT=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

public:
	// Required. The unique ID of the FORMAT field. Examples include "GT", "PL".
	std::string id;

	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The number of entries expected. See description above in the
	// VcfInfo message.
	std::string number;

	// Required. The type of the field. Valid values are "Integer", "Float",
	// "Character", and "String" (same as INFO except "Flag" is not supported).
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;
};

struct VcfFilter{
public:
	VcfFilter() : idx(0){}
	~VcfFilter() = default;

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FILTER=<ID=" + this->id;
		ret += ",Description=" + this->description;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FILTER=<ID=" + this->id;
		ret += ",Description=" + this->description;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

	friend std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt){
		stream.write((const char*)&flt.idx, sizeof(uint32_t));

		utility::SerializeString(flt.id, stream);
		utility::SerializeString(flt.description, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfFilter& flt){
		stream.read((char*)&flt.idx, sizeof(uint32_t));

		utility::DeserializeString(flt.id, stream);
		utility::DeserializeString(flt.description, stream);

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfFilter& flt){
		io::SerializePrimitive(flt.idx, buffer);
		io::SerializeString(flt.id, buffer);
		io::SerializeString(flt.description, buffer);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfFilter& flt){
		io::DeserializePrimitive(flt.idx, buffer);
		io::DeserializeString(flt.id, buffer);
		io::DeserializeString(flt.description, buffer);
		return(buffer);
	}

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the filter. Examples include "PASS", "RefCall".
	std::string id;

	// Required by VCF. The description of the filter.
	std::string description;
};

// This record type is a catch-all for other types of headers. For example,
// ##pedigreeDB=http://url_of_pedigrees
// The VcfExtra message would represent this with key="pedigreeDB",
// value="http://url_of_pedigrees".
struct VcfExtra{
public:
	VcfExtra() = default;
	VcfExtra(const std::string& key, const std::string& value) :
		key(key),
		value(value)
	{}

	~VcfExtra() = default;

	std::string ToVcfString(void) const{
		// Template:
		// ##source=CombineGVCFs
		std::string ret = "##" + this->key + "=" + this->value;
		return(ret);
	}

	friend std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra){
		utility::SerializeString(extra.key, stream);
		utility::SerializeString(extra.value, stream);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfExtra& extra){
		utility::DeserializeString(extra.key, stream);
		utility::DeserializeString(extra.value, stream);
		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfExtra& extra){
		io::SerializeString(extra.key, buffer);
		io::SerializeString(extra.value, buffer);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfExtra& extra){
		io::DeserializeString(extra.key, buffer);
		io::DeserializeString(extra.value, buffer);
		return(buffer);
	}

public:
  // Required by VCF. The key of the extra header field. Note that this key does
  // not have to be unique within a VcfHeader.
  std::string key;

  // Required by VCF. The value of the extra header field.
  std::string value;
};

// This record type is a catch-all for other headers containing multiple
// key-value pairs. For example, headers may have META lines that provide
// metadata about the VCF as a whole, e.g.
// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
// The VcfStructuredExtra message would represent this with key="META",
// and fields mapping "ID" -> "Assay", "Type" -> "String", etc.
struct VcfStructuredExtra{
public:
	VcfStructuredExtra() = default;
	~VcfStructuredExtra() = default;

	std::string ToVcfString(void) const{
		// Template:
		// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
		std::string ret = "##" + this->key + "=<";
		ret += this->fields[0].key + "=" + this->fields[0].value;
		for(uint32_t i = 1; i < this->fields.size(); ++i)
			ret += "," + this->fields[i].key + "=" + this->fields[i].value;
		ret += ">";
		return(ret);
	}

	friend std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra){
		utility::SerializeString(extra.key, stream);
		size_t l_extra = extra.fields.size();
		stream.write((const char*)&l_extra, sizeof(size_t));
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			stream << extra.fields[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra){
		utility::DeserializeString(extra.key, stream);
		size_t l_extra;
		stream.read((char*)&l_extra, sizeof(size_t));
		extra.fields.resize(l_extra);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			stream >> extra.fields[i];

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfStructuredExtra& extra){
		io::SerializeString(extra.key, buffer);
		size_t l_extra = extra.fields.size();
		io::SerializePrimitive(l_extra, buffer);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			buffer << extra.fields[i];

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfStructuredExtra& extra){
		io::DeserializeString(extra.key, buffer);
		size_t l_extra;
		io::DeserializePrimitive(l_extra, buffer);
		extra.fields.resize(l_extra);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			buffer >> extra.fields[i];

		return(buffer);
	}

public:
	// Required by VCF. The key of the extra header field. Note that this key does
	// not have to be unique within a VcfHeader.
	std::string key;

	// Required by VCF. The key=value pairs contained in the structure.
	std::vector<VcfExtra> fields;
};

struct YonContig : public VcfContig {
public:
	YonContig() : n_blocks(0){}
	YonContig(const VcfContig& vcf_contig) : VcfContig(vcf_contig), n_blocks(0){}
	~YonContig() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(VcfContig::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(VcfContig::ToVcfString(idx)); }

	inline void operator++(void){ ++this->n_blocks; }
	inline void operator--(void){ --this->n_blocks; }
	template <class T> inline void operator+=(const T value){ this->n_blocks += value; }
	template <class T> inline void operator-=(const T value){ this->n_blocks -= value; }

	friend std::ostream& operator<<(std::ostream& stream, const YonContig& contig){
		utility::SerializePrimitive(contig.idx, stream);
		utility::SerializePrimitive(contig.n_bases, stream);
		utility::SerializePrimitive(contig.n_blocks, stream);
		utility::SerializeString(contig.name, stream);
		utility::SerializeString(contig.description, stream);

		size_t size_helper = contig.extra.size();
		utility::SerializePrimitive(size_helper, stream);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			utility::SerializeString(contig.extra[i].first, stream);
			utility::SerializeString(contig.extra[i].second, stream);
		}
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonContig& contig){
		utility::DeserializePrimitive(contig.idx, stream);
		utility::DeserializePrimitive(contig.n_bases, stream);
		utility::DeserializePrimitive(contig.n_blocks, stream);
		utility::DeserializeString(contig.name, stream);
		utility::DeserializeString(contig.description, stream);

		size_t l_extra;
		utility::DeserializePrimitive(l_extra, stream);
		contig.extra.resize(l_extra);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			utility::DeserializeString(contig.extra[i].first, stream);
			utility::DeserializeString(contig.extra[i].second, stream);
		}
		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonContig& contig){
		io::SerializePrimitive(contig.idx, buffer);
		io::SerializePrimitive(contig.n_bases, buffer);
		io::SerializePrimitive(contig.n_blocks, buffer);
		io::SerializeString(contig.name, buffer);
		io::SerializeString(contig.description, buffer);

		size_t size_helper = contig.extra.size();
		io::SerializePrimitive(size_helper, buffer);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			io::SerializeString(contig.extra[i].first, buffer);
			io::SerializeString(contig.extra[i].second, buffer);
		}
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonContig& contig){
		io::DeserializePrimitive(contig.idx, buffer);
		io::DeserializePrimitive(contig.n_bases, buffer);
		io::DeserializePrimitive(contig.n_blocks, buffer);
		io::DeserializeString(contig.name, buffer);
		io::DeserializeString(contig.description, buffer);

		size_t l_extra;
		io::DeserializePrimitive(l_extra, buffer);
		contig.extra.resize(l_extra);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			io::DeserializeString(contig.extra[i].first, buffer);
			io::DeserializeString(contig.extra[i].second, buffer);
		}
		return(buffer);
	}

public:
	// Number of Tachyon blocks associated with this contig
	uint32_t n_blocks;
};

class YonInfo : public VcfInfo {
public:
	YonInfo() : yon_type(YON_VCF_HEADER_FLAG){}
	YonInfo(const VcfInfo& vcf_info) : VcfInfo(vcf_info), yon_type(YON_VCF_HEADER_FLAG){
		this->EvaluateType();
	}
	~YonInfo() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(VcfInfo::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(VcfInfo::ToVcfString(idx)); }

	bool EvaluateType(void){
		if(this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
		else if(this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
		else if(this->type == "Flag") this->yon_type = YON_VCF_HEADER_FLAG;
		else if(this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
		else if(this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
		else {
			std::cerr << "Illegal header type: " << this->type << std::endl;
			return false;
 		}
 		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const YonInfo& info){
		utility::SerializePrimitive(info.idx, stream);
		utility::SerializeString(info.id, stream);
		utility::SerializeString(info.number, stream);
		utility::SerializeString(info.type, stream);
		utility::SerializeString(info.description, stream);
		utility::SerializeString(info.source, stream);
		utility::SerializeString(info.version, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonInfo& info){
		utility::DeserializePrimitive(info.idx, stream);
		utility::DeserializeString(info.id, stream);
		utility::DeserializeString(info.number, stream);
		utility::DeserializeString(info.type, stream);
		utility::DeserializeString(info.description, stream);
		utility::DeserializeString(info.source, stream);
		utility::DeserializeString(info.version, stream);
		info.EvaluateType();

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonInfo& info){
		io::SerializePrimitive(info.idx, buffer);
		io::SerializeString(info.id, buffer);
		io::SerializeString(info.number, buffer);
		io::SerializeString(info.type, buffer);
		io::SerializeString(info.description, buffer);
		io::SerializeString(info.source, buffer);
		io::SerializeString(info.version, buffer);

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonInfo& info){
		io::DeserializePrimitive(info.idx, buffer);
		io::DeserializeString(info.id, buffer);
		io::DeserializeString(info.number, buffer);
		io::DeserializeString(info.type, buffer);
		io::DeserializeString(info.description, buffer);
		io::DeserializeString(info.source, buffer);
		io::DeserializeString(info.version, buffer);
		info.EvaluateType();

		return(buffer);
	}

public:
	TACHYON_VARIANT_HEADER_FIELD_TYPE yon_type;
};

class YonFormat : public VcfFormat {
public:
	YonFormat() : yon_type(YON_VCF_HEADER_FLAG){}
	YonFormat(const VcfFormat& vcf_format) : VcfFormat(vcf_format), yon_type(YON_VCF_HEADER_FLAG){
		this->EvaluateType();
	}
	~YonFormat() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(VcfFormat::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(VcfFormat::ToVcfString(idx)); }

	bool EvaluateType(void){
		if(this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
		else if(this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
		else if(this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
		else if(this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
		else {
			std::cerr << "Illegal header type: " << this->type << std::endl;
			return false;
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const YonFormat& fmt){
		utility::SerializePrimitive(fmt.idx, stream);
		utility::SerializeString(fmt.id, stream);
		utility::SerializeString(fmt.number, stream);
		utility::SerializeString(fmt.type, stream);
		utility::SerializeString(fmt.description, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonFormat& fmt){
		utility::DeserializePrimitive(fmt.idx, stream);
		utility::DeserializeString(fmt.id, stream);
		utility::DeserializeString(fmt.number, stream);
		utility::DeserializeString(fmt.type, stream);
		utility::DeserializeString(fmt.description, stream);
		fmt.EvaluateType();

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonFormat& fmt){
		io::SerializePrimitive(fmt.idx, buffer);
		io::SerializeString(fmt.id, buffer);
		io::SerializeString(fmt.number, buffer);
		io::SerializeString(fmt.type, buffer);
		io::SerializeString(fmt.description, buffer);

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonFormat& fmt){
		io::DeserializePrimitive(fmt.idx, buffer);
		io::DeserializeString(fmt.id, buffer);
		io::DeserializeString(fmt.number, buffer);
		io::DeserializeString(fmt.type, buffer);
		io::DeserializeString(fmt.description, buffer);
		fmt.EvaluateType();

		return(buffer);
	}

public:
	TACHYON_VARIANT_HEADER_FIELD_TYPE yon_type;
};

namespace utility {

/**<
 * Helper functions for converting a given encoded tachyon primitive value
 * into a valid Vcf string.
 * @param stream Dst buffer reference.
 * @param data   Src pointer to data.
 * @param n_data Number of entries the src pointer is pointing to.
 * @return       Returns a reference to the dst buffer reference.
 */
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const uint8_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const uint16_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const uint32_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const uint64_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const int8_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const int16_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const int32_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const int64_t* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const char* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const float* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const double* const data, const size_t n_data);
io::BasicBuffer& ToVcfString(io::BasicBuffer& stream, const std::string& string);

/**<
 * Helper functions for adding raw primitive data into a htslib bcf1_t info field.
 * @param rec       Dst htslib bcf1_t record pointer.
 * @param hdr       Pointer reference to htslib bcf header.
 * @param tag       Tag string telling us where to append data in the record.
 * @param data      Src data pointer.
 * @param n_entries Number of entries the src data pointer is pointing to.
 * @return          Returns the (possibly) updated htslib bcf1_t record pointer.
 */
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint8_t* const data,  const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint16_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint32_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint64_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int8_t* const data,   const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int16_t* const data,  const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int32_t* const data,  const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const float* const data,    const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const double* const data,   const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const std::string& data);

/**<
 * Conversion functions for a given primitive into the accepted primitive types
 * of htslib bcf1_t. Pre-processing steps involves the narrowing or expansion of
 * integers/floats and their corresponding missing/NAN/EOV values.
 * @param src       Src data pointer.
 * @param dst       Dst primitive pointer.
 * @param n_entries Number of entries the src data pointer is pointing to.
 * @return          Returns the (possibly) converted primitive.
 */
int32_t* FormatDataHtslib(const uint8_t*  const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint16_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint32_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint64_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int8_t*  const src,  int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int16_t* const src,  int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int32_t* const src,  int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int64_t* const src,  int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const float*   const src,  int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const double*  const src,  int32_t* dst, const size_t n_entries);

template <class T>
float* FormatDataHtslib(const T* const src, float* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}


}
}

#endif /* UTILITY_SUPPORT_VCF_H_ */
