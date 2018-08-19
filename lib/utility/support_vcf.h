#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_


#include <iostream>
#include <cmath>

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
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const uint8_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const uint16_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const uint32_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const uint64_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const int8_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const int16_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const int32_t* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const char* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const float* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const double* const data, const size_t n_data);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const std::string& string);

std::ostream& to_vcf_string(std::ostream& stream, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);

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
