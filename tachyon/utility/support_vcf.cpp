#include "support_vcf.h"

namespace tachyon{
namespace util{

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<BYTE>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	stream << container[0];
	for(U32 i = 1; i < container.size(); ++i)
		stream << ',' << container[i];

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U16>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	stream << container[0];
	for(U32 i = 1; i < container.size(); ++i)
		stream << ',' << container[i];

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U32>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	stream << container[0];
	for(U32 i = 1; i < container.size(); ++i)
		stream << ',' << container[i];

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<U64>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	stream << container[0];
	for(U32 i = 1; i < container.size(); ++i)
		stream << ',' << container[i];

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<SBYTE>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	const BYTE* const ref = reinterpret_cast<const BYTE* const>(container.data());

	// If the first value is end-of-vector then return
	if(ref[0] == YON_BYTE_EOV)
		return(stream.put('.'));

	// First value
	if(ref[0] == YON_BYTE_MISSING) stream << '.';
	else stream << container[0];

	// Remainder values
	for(U32 i = 1; i < container.size(); ++i){
		if(ref[i] == YON_BYTE_MISSING) stream << ",.";
		else if(ref[i] == YON_BYTE_EOV){ return stream; }
		else stream << ',' << container[i];
	}

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<S16>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	const U16* const ref = reinterpret_cast<const U16* const>(container.data());

	// If the first value is end-of-vector then return
	if(ref[0] == YON_SHORT_EOV)
		return(stream.put('.'));

	// First value
	if(ref[0] == YON_SHORT_MISSING) stream << '.';
	else stream << container[0];

	// Remainder values
	for(U32 i = 1; i < container.size(); ++i){
		if(ref[i] == YON_SHORT_MISSING) stream << ",.";
		else if(ref[i] == YON_SHORT_EOV){ return stream; }
		else stream << ',' << container[i];
	}

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<S32>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	const U32* const ref = reinterpret_cast<const U32* const>(container.data());

	// If the first value is end-of-vector then return
	if(ref[0] == YON_INT_EOV)
		return(stream.put('.'));

	// First value
	if(ref[0] == YON_INT_MISSING) stream << '.';
	else stream << container[0];

	// Remainder values
	for(U32 i = 1; i < container.size(); ++i){
		if(ref[i] == YON_INT_MISSING) stream << ",.";
		else if(ref[i] == YON_INT_EOV){ return stream; }
		else stream << ',' << container[i];
	}

	return(stream);
}

// Special case
std::ostream& to_vcf_string_char(std::ostream& stream, containers::PrimitiveContainer<char>& container){
	if(container.size() == 0)
		return(stream.put('.'));

	stream << container[0];
	for(U32 i = 1; i < container.size(); ++i)
		stream << ',' << container[i];

	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, containers::PrimitiveContainer<float>& container){
	if(container.size() == 0)
		return(stream.put('.'));

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


}
}
