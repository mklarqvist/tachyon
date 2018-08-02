#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_


#include <iostream>
#include <cmath>
#include "support/type_definitions.h"
//#include "containers/primitive_container.h"
#include "core/genotype_object.h"
#include "core/data_block_settings.h"

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

// Genotype objects
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const core::GTObject& gt_object);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& stream, const std::vector<core::GTObject>& gt_objects);

std::ostream& to_vcf_string(std::ostream& stream, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header);
io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);
io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller);

}
}



#endif /* UTILITY_SUPPORT_VCF_H_ */
