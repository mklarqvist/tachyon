#include "support_vcf.h"

namespace tachyon {
namespace utility {

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint8_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble((uint32_t)data[0]);
	for(uint32_t i = 1; i < n_data; ++i){
		buffer += ',';
		buffer.AddReadble((uint32_t)data[i]);
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint16_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for(uint32_t i = 1; i < n_data; ++i){
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint32_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for(uint32_t i = 1; i < n_data; ++i){
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint64_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for(uint32_t i = 1; i < n_data; ++i){
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const int8_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if(data[0] == YON_BYTE_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(data[0] == YON_BYTE_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(data[i] == YON_BYTE_MISSING) buffer += ",.";
		else if(data[i] == YON_BYTE_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const int16_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}
	// If the first value is end-of-vector then return
	if(data[0] == YON_SHORT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(data[0] == YON_SHORT_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(data[i] == YON_SHORT_MISSING) buffer += ",.";
		else if(data[i] == YON_SHORT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const int32_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if(data[0] == YON_INT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(data[0] == YON_INT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(data[i] == YON_INT_MISSING) buffer += ",.";
		else if(data[i] == YON_INT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const int64_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if(data[0] == YON_INT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(data[0] == YON_INT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(data[i] == YON_INT_MISSING) buffer += ",.";
		else if(data[i] == YON_INT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

// Special case
io::BasicBuffer& ToVcfString_char(io::BasicBuffer& buffer, const char* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	buffer += data[0];
	for(uint32_t i = 1; i < n_data; ++i){
		buffer += ',';
		buffer += data[i];
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const float* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	const uint32_t* const ref = reinterpret_cast<const uint32_t* const>(data);

	// If the first value is end-of-vector then return
	if(ref[0] == YON_FLOAT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(ref[0] == YON_FLOAT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(ref[i] == YON_FLOAT_MISSING) buffer += ",.";
		else if(ref[i] == YON_FLOAT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const double* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	const uint32_t* const ref = reinterpret_cast<const uint32_t* const>(data);

	// If the first value is end-of-vector then return
	if(ref[0] == YON_FLOAT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(ref[0] == YON_FLOAT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(ref[i] == YON_FLOAT_MISSING) buffer += ",.";
		else if(ref[i] == YON_FLOAT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

int32_t* FormatDataHtslib(const uint8_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint16_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint32_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint64_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const int8_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i){
		if(src[i] == INT8_MIN)  { dst[i] = bcf_int32_missing; }
		else if(src[i] == INT8_MIN+1){ dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int16_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i){
		if(src[i] == INT16_MIN)  { dst[i] = bcf_int32_missing; }
		else if(src[i] == INT16_MIN+1){ dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int32_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i){
		if(src[i] == INT32_MIN)  { dst[i] = bcf_int32_missing; }
		else if(src[i] == INT32_MIN+1){ dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int64_t* const src, int32_t* dst, const size_t n_entries){
	for(uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const float* const src, int32_t* dst, const size_t n_entries){
	std::cerr << utility::timestamp("ERROR") << "Cannot convert float into integer." << std::endl;
	return(dst);
}

int32_t* FormatDataHtslib(const double* const src, int32_t* dst, const size_t n_entries){
	std::cerr << utility::timestamp("ERROR") << "Cannot convert double into integer." << std::endl;
	return(dst);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const uint8_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), &tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const uint16_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const uint32_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const uint64_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const int8_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const int16_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const int32_t* const data,
                                  const size_t n_entries)
{
	bcf_update_info_int32(hdr, rec, tag.data(), data, n_entries);
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const int64_t* const data,
                                  const size_t n_entries)
{
	int32_t* tmpi = new int32_t[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_int32(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const float* const data,
                                  const size_t n_entries)
{
	float* tmpi = new float[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_float(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const double* const data,
                                  const size_t n_entries)
{
	float* tmpi = new float[n_entries];
	FormatDataHtslib(data, tmpi, n_entries);
	bcf_update_info_float(hdr, rec, tag.data(), tmpi, n_entries);
	delete [] tmpi;
	return(rec);
}

bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                  bcf_hdr_t* hdr,
                                  const std::string& tag,
                                  const std::string& data)
{
	bcf_update_info_string(hdr, rec, tag.data(), data.data());
	return(rec);
}

}
}
