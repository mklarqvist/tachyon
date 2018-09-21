#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_

#include <iostream>
#include <cmath>

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

#include "core/data_block_settings.h"

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
struct VcfInfo {
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
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfFilter& flt);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfFilter& flt);

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
struct VcfExtra {
public:
	VcfExtra() = default;
	VcfExtra(const std::string& key, const std::string& value);
	~VcfExtra() = default;

	std::string ToVcfString(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra);
	friend std::istream& operator>>(std::istream& stream, VcfExtra& extra);
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfExtra& extra);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfExtra& extra);

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
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfStructuredExtra& extra);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfStructuredExtra& extra);

public:
	// Required by VCF. The key of the extra header field. Note that this key does
	// not have to be unique within a VcfHeader.
	std::string key;

	// Required by VCF. The key=value pairs contained in the structure.
	std::vector<VcfExtra> fields;
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
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint8_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint16_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint32_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const uint64_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int8_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int16_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const int32_t* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const float* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const double* const data, const size_t n_entries);
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
int32_t* FormatDataHtslib(const uint8_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint16_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint32_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const uint64_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int8_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int16_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int32_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const int64_t* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const float* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const double* const src, int32_t* dst, const size_t n_entries);

template <class T>
float* FormatDataHtslib(const T* const src, float* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}


}
}

#endif /* UTILITY_SUPPORT_VCF_H_ */
