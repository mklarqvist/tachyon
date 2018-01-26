#ifndef UTILITY_SUPPORT_VCF_H_
#define UTILITY_SUPPORT_VCF_H_


#include <iostream>
#include <cmath>
#include "../support/TypeDefinitions.h"
#include "../containers/primitive_container.h"

namespace tachyon{
namespace util{

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

// Primtive container declarations
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<BYTE>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<U16>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<U32>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<U64>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<SBYTE>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<S16>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<S32>& container);
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<char>& container);
inline std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<float>& container){
	if(container.size() == 0)
		return(stream);

	const U32* const ref = reinterpret_cast<const U32* const>(container.data());

	// If the first value is end-of-vector then return
	if(ref[0] == YON_FLOAT_EOV)
		return(stream.put('.'));

	// First value
	if(ref[0] == YON_FLOAT_MISSING) stream << '.';
	else stream << container[0];

	// Remainder values
	for(U32 i = 1; i < container.size(); ++i){
		if(ref[i] == YON_FLOAT_MISSING) stream << ",.";
		else if(ref[i] == YON_FLOAT_EOV){ return stream; }
		else stream << ',' << container[i];
	}

	return(stream);
}
std::ostream& to_vcf_string(std::ostream& stream, core::PrimitiveContainer<double>& container);

}
}



#endif /* UTILITY_SUPPORT_VCF_H_ */
