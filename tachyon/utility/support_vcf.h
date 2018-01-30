#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_


#include <iostream>
#include <cmath>
#include "../support/type_definitions.h"
#include "../containers/primitive_container.h"

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
std::ostream& to_vcf_string(std::ostream& stream, const BYTE* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const U16* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const U32* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const U64* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const SBYTE* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const S16* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const S32* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const char* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const float* const values, const U32 n_entries);
std::ostream& to_vcf_string(std::ostream& stream, const double* const values, const U32 n_entries);

// Primitive container declarations
// Unsigned values does not have missing values
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<BYTE>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U16>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U32>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U64>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<SBYTE>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<S16>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<S32>& container);
// Special case
std::ostream& to_vcf_string_char(std::ostream& stream, containers::PrimitiveContainer<char>& container);
std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<float>& container);

//std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<double>& container);

}
}



#endif /* UTILITY_SUPPORT_VCF_H_ */
