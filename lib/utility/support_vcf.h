#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_


#include <iostream>
#include <cmath>
#include "support/type_definitions.h"
//#include "containers/primitive_container.h"
#include "core/data_block_settings.h"
#include "core/meta_entry.h"

namespace tachyon{
namespace utility{

#define YON_BYTE_MISSING        0x80
#define YON_BYTE_EOV            0x81
#define YON_SHORT_MISSING       0x8000
#define YON_SHORT_EOV           0x8001
#define YON_INT_MISSING         0x80000000
#define YON_INT_EOV             0x80000001
#define YON_FLOAT_NAN           0x7FC00000
#define YON_FLOAT_MISSING       0x7F800001
#define YON_FLOAT_EOV           0x7F800002

// Base functionality converting data to a valid VCF string
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const BYTE* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const U16* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const U32* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const U64* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const SBYTE* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const S16* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const S32* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const char* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const float* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const double* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const std::string& string);

std::ostream& to_vcf_string(std::ostream& stream, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const BYTE* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const U16* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const U32* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const U64* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const SBYTE* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const S16* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const S32* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const float* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const double* const data, const size_t n_entries);
bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag, const std::string& data);

int32_t* FormatDataHtslib(const BYTE* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const U16* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const U32* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const U64* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const SBYTE* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const S16* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const S32* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const S64* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const float* const src, int32_t* dst, const size_t n_entries);
int32_t* FormatDataHtslib(const double* const src, int32_t* dst, const size_t n_entries);

template <class T>
float* FormatDataHtslib(const T* const src, float* dst, const size_t n_entries){
	for(U32 i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}


}
}

#endif /* UTILITY_SUPPORT_VCF_H_ */
