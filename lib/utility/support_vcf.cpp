#include "support_vcf.h"

namespace tachyon{
namespace utility{

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const uint8_t* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const uint16_t* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const uint32_t* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const uint64_t* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const int8_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	const uint8_t* const ref = reinterpret_cast<const uint8_t* const>(data);

	// If the first value is end-of-vector then return
	if(ref[0] == YON_BYTE_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(ref[0] == YON_BYTE_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(ref[i] == YON_BYTE_MISSING) buffer += ",.";
		else if(ref[i] == YON_BYTE_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const int16_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	const uint16_t* const ref = reinterpret_cast<const uint16_t* const>(data);

	// If the first value is end-of-vector then return
	if(ref[0] == YON_SHORT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(ref[0] == YON_SHORT_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(ref[i] == YON_SHORT_MISSING) buffer += ",.";
		else if(ref[i] == YON_SHORT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const int32_t* data, const size_t n_data){
	if(n_data == 0){
		buffer += '.';
		return(buffer);
	}

	const uint32_t* const ref = reinterpret_cast<const uint32_t* const>(data);

	// If the first value is end-of-vector then return
	if(ref[0] == YON_INT_EOV){
		buffer += '.';
		return(buffer);
	}

	// First value
	if(ref[0] == YON_INT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for(uint32_t i = 1; i < n_data; ++i){
		if(ref[i] == YON_INT_MISSING) buffer += ",.";
		else if(ref[i] == YON_INT_EOV){ return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

// Special case
io::BasicBuffer& to_vcf_string_char(io::BasicBuffer& buffer, const char* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const float* data, const size_t n_data){
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

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const double* data, const size_t n_data){
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

std::ostream& to_vcf_string(std::ostream& stream, const std::string& string){
	stream << string;
	return(stream);
}

std::ostream& to_vcf_string(std::ostream& stream, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header){
	stream.write(&header.GetContig(meta_entry.contigID)->name[0], header.GetContig(meta_entry.contigID)->name.size()) << '\t';
	stream << meta_entry.position + 1 << delimiter;

	if(meta_entry.name.size() == 0) stream.put('.');
	else stream.write(&meta_entry.name[0], meta_entry.name.size());
	stream.put(delimiter);
	if(meta_entry.n_alleles){
		stream.write(meta_entry.alleles[0].allele, meta_entry.alleles[0].l_allele);
		//stream << meta_entry.alleles[0].l_allele;
		stream.put(delimiter);
		stream.write(meta_entry.alleles[1].allele, meta_entry.alleles[1].l_allele);
		for(uint32_t i = 2; i < meta_entry.n_alleles; ++i){
			stream.put(',');
			stream.write(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
		}
	} else {
		//stream << ".\t.\t";
		stream << '.' << delimiter << '.' << delimiter;
	}

	if(std::isnan(meta_entry.quality)){
		stream << delimiter << '.' << delimiter;
		//stream << "\t.\t";
	}
	else {
		stream << delimiter << meta_entry.quality << delimiter;
		//stream << '\t' << meta_entry.quality << '\t';
	}
	return(stream);
}

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header){
	buffer += header.GetContig(meta_entry.contigID)->name;
	buffer += delimiter;
	buffer.AddReadble(meta_entry.position + 1);
	buffer += delimiter;

	if(meta_entry.name.size() == 0) buffer += '.';
	else buffer += meta_entry.name;
	buffer += delimiter;
	if(meta_entry.n_alleles){
		buffer.Add(meta_entry.alleles[0].allele, meta_entry.alleles[0].l_allele);
		buffer += delimiter;
		buffer.Add(meta_entry.alleles[1].allele, meta_entry.alleles[1].l_allele);
		for(uint32_t i = 2; i < meta_entry.n_alleles; ++i){
			buffer += ',';
			buffer.Add(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
		}
		buffer += delimiter;
	} else {
		buffer += '.';
		buffer += delimiter;
		buffer += '.';
		buffer += delimiter;
	}

	if(std::isnan(meta_entry.quality)){
		buffer += '.';
	}
	else {
		buffer.AddReadble(meta_entry.quality);
	}

	return(buffer);
}

io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller){
	//if(controller.contig.display){
		buffer += header.GetContig(meta_entry.contigID)->name;
		buffer += delimiter;
	//}

	//if(controller.positions.display){
		buffer.AddReadble(meta_entry.position + 1);
		buffer += delimiter;
	//}

	//if(controller.names.display){
		if(meta_entry.name.size() == 0) buffer += '.';
		else buffer += meta_entry.name;
		buffer += delimiter;
	//}

	if(controller.display_ref){
		if(meta_entry.n_alleles) buffer.Add(meta_entry.alleles[0].allele, meta_entry.alleles[0].l_allele);
		else buffer += '.';
		buffer += delimiter;
	}

	if(controller.display_alt){
		if(meta_entry.n_alleles){
			buffer.Add(meta_entry.alleles[1].allele, meta_entry.alleles[1].l_allele);
			for(uint32_t i = 2; i < meta_entry.n_alleles; ++i){
				buffer += ',';
				buffer.Add(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
			}
		} else
			buffer += '.';

		buffer += delimiter;
	}

	//if(controller.quality.display){
		if(std::isnan(meta_entry.quality)) buffer += '.';
		else buffer.AddReadble(meta_entry.quality);
		buffer += delimiter;
	//}

	return(buffer);
}

io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller){
	return(utility::to_json_string(buffer, meta_entry, header, controller));
}


io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const char& delimiter, const core::MetaEntry& meta_entry, const VariantHeader& header, const DataBlockSettings& controller){
	bool add = false;
	//if(controller.contig.display){
		buffer += "\"contig\":\"";
		buffer += header.GetContig(meta_entry.contigID)->name;
		buffer += '"';
		add = true;
	//}

	//if(controller.positions.display){
		if(add){ buffer += ','; add = false; }
		buffer += "\"position\":";
		buffer.AddReadble(meta_entry.position + 1);
		add = true;
	//}

	//if(controller.names.display){
		if(add){ buffer += ','; add = false; }
		buffer += "\"name\":";
		if(meta_entry.name.size() == 0) buffer += "null";
		else {
			buffer += '"';
			buffer += meta_entry.name;
			buffer += '"';
			add = true;
		}
	//}

	//if(controller.display_ref){
		if(add){ buffer += ','; add = false; }
		buffer += "\"ref\":";
		if(meta_entry.n_alleles){
			buffer += '"';
			buffer.Add(meta_entry.alleles[0].allele, meta_entry.alleles[0].l_allele);
			buffer += '"';
		} else {
			buffer += "null";
		}
		add = true;
	//}

	if(controller.display_alt){
		if(add){ buffer += ','; add = false; }
		buffer += "\"alt\":";
		if(meta_entry.n_alleles){
			if(meta_entry.n_alleles > 2) buffer += '[';
			buffer += '"';
			buffer.Add(meta_entry.alleles[1].allele, meta_entry.alleles[1].l_allele);
			buffer += '"';
			for(uint32_t i = 2; i < meta_entry.n_alleles; ++i){
				buffer += ',';
				buffer += '"';
				buffer.Add(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
				buffer += '"';
			}
			if(meta_entry.n_alleles > 2) buffer += ']';
		} else {
			buffer += "null";
		}
		add = true;
	}

	//if(controller.quality.display){
		if(add){ buffer += ','; add = false; }
		buffer += "\"quality\":";
		if(std::isnan(meta_entry.quality)) buffer += 0;
		else buffer.AddReadble(meta_entry.quality);
		add = true;
	//}

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
