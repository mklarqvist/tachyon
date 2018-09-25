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
#ifndef TACHYON_SUPPORT_VCF_H_
#define TACHYON_SUPPORT_VCF_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "tachyon.h"
#include "buffer.h"

#include "htslib/vcf.h"

#define YON_BYTE_MISSING  INT8_MIN
#define YON_BYTE_EOV      (INT8_MIN+1)
#define YON_SHORT_MISSING INT16_MIN
#define YON_SHORT_EOV     (INT16_MIN+1)
#define YON_INT_MISSING   INT32_MIN
#define YON_INT_EOV       (INT32_MIN+1)
#define YON_FLOAT_NAN     0x7FC00000
#define YON_FLOAT_MISSING 0x7F800001
#define YON_FLOAT_EOV     0x7F800002
#define YON_BCF_GT_MISSING 0

// taken from htslib vcf.c and renamed for convenience
static inline int yon_float_is_missing(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==YON_FLOAT_MISSING ? 1 : 0;
}

static inline int yon_float_is_vector_end(float f)
{
    union { uint32_t i; float f; } u;
    u.f = f;
    return u.i==YON_FLOAT_EOV ? 1 : 0;
}

namespace tachyon {

//
// -----------------------------------------------------------------------------
// VCF type encoding utilities
template<class T>
struct VcfType {
  // Predicates for checking missing and sentinel entries.  Use these, not ==.
  // Is argument the "missing" value?
  static bool IsMissing(T);
  // Is argument the vector end sentinel value?
  static bool IsVectorEnd(T);
};

// See interface description comment above.
template<>
struct VcfType<int8_t> {
  static bool IsMissing(int8_t v)  { return (v == YON_BYTE_MISSING); }
  static bool IsVectorEnd(int8_t v){ return (v == YON_BYTE_EOV); }
};

// See interface description comment above.
template<>
struct VcfType<int16_t> {
  static bool IsMissing(int16_t v)  { return (v == YON_SHORT_MISSING); }
  static bool IsVectorEnd(int16_t v){ return (v == YON_SHORT_EOV); }
};

// See interface description comment above.
template<>
struct VcfType<int> {
  static bool IsMissing(int v)  { return (v == YON_INT_MISSING); }
  static bool IsVectorEnd(int v){ return (v == YON_INT_EOV); }
};

template <class T>
struct VcfGenotype {
	// Predicates for checking missing and sentinel entries.
	static bool IsMissing(const T& value){ return(value == YON_BCF_GT_MISSING); }
};

template<>
struct VcfType<float> {
  static bool IsMissing(float v)  { return yon_float_is_missing(v); }
  static bool IsVectorEnd(float v){ return yon_float_is_vector_end(v); }
};


struct VcfContig {
public:
	VcfContig();
	~VcfContig() = default;

	std::string ToVcfString(const bool is_bcf = false) const;

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
	VcfInfo();
	~VcfInfo() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

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
	VcfFormat();
	~VcfFormat() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

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
	VcfFilter();
	~VcfFilter() = default;

	std::string ToVcfString(const bool is_bcf = false) const;
	std::string ToVcfString(const uint32_t idx) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt);
	friend std::istream& operator>>(std::istream& stream, VcfFilter& flt);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfFilter& flt);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfFilter& flt);

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
	VcfExtra(const std::string& key, const std::string& value);
	~VcfExtra() = default;

	std::string ToVcfString(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra);
	friend std::istream& operator>>(std::istream& stream, VcfExtra& extra);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfExtra& extra);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfExtra& extra);

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

	std::string ToVcfString(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra);
	friend std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfStructuredExtra& extra);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfStructuredExtra& extra);

public:
	// Required by VCF. The key of the extra header field. Note that this key does
	// not have to be unique within a VcfHeader.
	std::string key;

	// Required by VCF. The key=value pairs contained in the structure.
	std::vector<VcfExtra> fields;
};

struct YonContig : public VcfContig {
public:
	YonContig();
	YonContig(const VcfContig& vcf_contig);
	~YonContig() = default;

	inline std::string ToVcfString(const bool is_bcf = false) const{ return(VcfContig::ToVcfString(is_bcf)); }
	inline std::string ToVcfString(const uint32_t idx) const{ return(VcfContig::ToVcfString(idx)); }

	inline void operator++(void){ ++this->n_blocks; }
	inline void operator--(void){ --this->n_blocks; }
	template <class T> inline void operator+=(const T value){ this->n_blocks += value; }
	template <class T> inline void operator-=(const T value){ this->n_blocks -= value; }

	friend std::ostream& operator<<(std::ostream& stream, const YonContig& contig);
	friend std::istream& operator>>(std::istream& stream, YonContig& contig);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonContig& contig);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, YonContig& contig);

public:
	// Number of Tachyon blocks associated with this contig
	uint32_t n_blocks;
};

class YonInfo : public VcfInfo {
public:
	YonInfo();
	YonInfo(const VcfInfo& vcf_info);
	~YonInfo() = default;

	inline std::string ToVcfString(const bool is_bcf = false) const{ return(VcfInfo::ToVcfString(is_bcf)); }
	inline std::string ToVcfString(const uint32_t idx) const{ return(VcfInfo::ToVcfString(idx)); }

	bool EvaluateType(void);

	friend std::ostream& operator<<(std::ostream& stream, const YonInfo& info);
	friend std::istream& operator>>(std::istream& stream, YonInfo& info);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonInfo& info);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, YonInfo& info);

public:
	TACHYON_VARIANT_HEADER_FIELD_TYPE yon_type;
};

class YonFormat : public VcfFormat {
public:
	YonFormat();
	YonFormat(const VcfFormat& vcf_format);
	~YonFormat() = default;

	inline std::string ToVcfString(const bool is_bcf = false) const{ return(VcfFormat::ToVcfString(is_bcf)); }
	inline std::string ToVcfString(const uint32_t idx) const{ return(VcfFormat::ToVcfString(idx)); }

	bool EvaluateType(void);

	friend std::ostream& operator<<(std::ostream& stream, const YonFormat& fmt);
	friend std::istream& operator>>(std::istream& stream, YonFormat& fmt);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonFormat& fmt);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, YonFormat& fmt);

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
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const uint8_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const uint16_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const uint32_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const uint64_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const int8_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const int16_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const int32_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const int64_t* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const char* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const float* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const double* const data, const size_t n_data);
yon_buffer_t& ToVcfString(yon_buffer_t& stream, const std::string& string);

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
