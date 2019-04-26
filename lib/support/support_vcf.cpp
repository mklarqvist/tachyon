#include "support_vcf.h"

namespace tachyon {

VcfContig::VcfContig() : idx(0), n_bases(0) {}

std::string VcfContig::ToVcfString(const bool is_bcf) const {
	// Template:
	// ##contig=<ID=GL000241.1,assembly=b37,length=42152>
	std::string ret = "##contig=<ID=" + this->name;
	if (extra.size()) {
		ret += "," + this->extra[0].first + "=" + this->extra[0].second;
		for (uint32_t i = 1; i < this->extra.size(); ++i) {
			ret += "," + this->extra[i].first + "=" + this->extra[i].second;
		}
	}
	if (this->description.size()) ret += ",Description=" + this->description;
	ret += ",length=" + std::to_string(this->n_bases);
	if (is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

VcfInfo::VcfInfo() : idx(0) {}

std::string VcfInfo::ToVcfString(const bool is_bcf) const {
	// Template:
	// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
	std::string ret = "##INFO=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if (this->source.size()) ret += ",Source=" + this->source;
	if (this->source.size()) ret += ",Version=" + this->version;
	if (is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfInfo::ToVcfString(const uint32_t idx) const {
	// Template:
	// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
	std::string ret = "##INFO=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if (this->source.size()) ret += ",Source=" + this->source;
	if (this->source.size()) ret += ",Version=" + this->version;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}

VcfFormat::VcfFormat() : idx(0) {}

std::string VcfFormat::ToVcfString(const bool is_bcf) const {
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FORMAT=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	if (is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfFormat::ToVcfString(const uint32_t idx) const {
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FORMAT=<ID=" + this->id;
	ret += ",Number=" + this->number;
	ret += ",Type=" + this->type;
	ret += ",Description=" + this->description;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}


VcfFilter::VcfFilter() : idx(0) {}

std::string VcfFilter::ToVcfString(const bool is_bcf) const {
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FILTER=<ID=" + this->id;
	ret += ",Description=" + this->description;
	if (is_bcf) ret += ",IDX=" + std::to_string(this->idx);
	ret += ">";
	return(ret);
}

std::string VcfFilter::ToVcfString(const uint32_t idx) const {
	// Template:
	// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
	std::string ret = "##FILTER=<ID=" + this->id;
	ret += ",Description=" + this->description;
	ret += ",IDX=" + std::to_string(idx);
	ret += ">";
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt) {
	stream.write((const char*)&flt.idx, sizeof(uint32_t));

	utility::SerializeString(flt.id, stream);
	utility::SerializeString(flt.description, stream);

	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfFilter& flt) {
	stream.read((char*)&flt.idx, sizeof(uint32_t));

	utility::DeserializeString(flt.id, stream);
	utility::DeserializeString(flt.description, stream);

	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfFilter& flt) {
	SerializePrimitive(flt.idx, buffer);
	SerializeString(flt.id, buffer);
	SerializeString(flt.description, buffer);
	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfFilter& flt) {
	DeserializePrimitive(flt.idx, buffer);
	DeserializeString(flt.id, buffer);
	DeserializeString(flt.description, buffer);
	return(buffer);
}

VcfExtra::VcfExtra(const std::string& key, const std::string& value) :
	key(key),
	value(value)
{}

std::string VcfExtra::ToVcfString(void) const {
	// Template:
	// ##source=CombineGVCFs
	std::string ret = "##" + this->key + "=" + this->value;
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra) {
	utility::SerializeString(extra.key, stream);
	utility::SerializeString(extra.value, stream);
	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfExtra& extra) {
	utility::DeserializeString(extra.key, stream);
	utility::DeserializeString(extra.value, stream);
	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfExtra& extra) {
	SerializeString(extra.key, buffer);
	SerializeString(extra.value, buffer);
	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfExtra& extra) {
	DeserializeString(extra.key, buffer);
	DeserializeString(extra.value, buffer);
	return(buffer);
}

std::string VcfStructuredExtra::ToVcfString(void) const {
	// Template:
	// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
	std::string ret = "##" + this->key + "=<";
	ret += this->fields[0].key + "=" + this->fields[0].value;
	for (uint32_t i = 1; i < this->fields.size(); ++i)
		ret += "," + this->fields[i].key + "=" + this->fields[i].value;
	ret += ">";
	return(ret);
}

std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra) {
	utility::SerializeString(extra.key, stream);
	size_t l_extra = extra.fields.size();
	stream.write((const char*)&l_extra, sizeof(size_t));
	for (uint32_t i = 0; i < extra.fields.size(); ++i)
		stream << extra.fields[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra) {
	utility::DeserializeString(extra.key, stream);
	size_t l_extra;
	stream.read((char*)&l_extra, sizeof(size_t));
	extra.fields.resize(l_extra);
	for (uint32_t i = 0; i < extra.fields.size(); ++i)
		stream >> extra.fields[i];

	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const VcfStructuredExtra& extra) {
	SerializeString(extra.key, buffer);
	size_t l_extra = extra.fields.size();
	SerializePrimitive(l_extra, buffer);
	for (uint32_t i = 0; i < extra.fields.size(); ++i)
		buffer << extra.fields[i];

	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, VcfStructuredExtra& extra) {
	DeserializeString(extra.key, buffer);
	size_t l_extra;
	DeserializePrimitive(l_extra, buffer);
	extra.fields.resize(l_extra);
	for (uint32_t i = 0; i < extra.fields.size(); ++i)
		buffer >> extra.fields[i];

	return(buffer);
}

YonContig::YonContig() : n_blocks(0) {}
YonContig::YonContig(const VcfContig& vcf_contig) : VcfContig(vcf_contig), n_blocks(0) {}

std::ostream& operator<<(std::ostream& stream, const YonContig& contig) {
	utility::SerializePrimitive(contig.idx, stream);
	utility::SerializePrimitive(contig.n_bases, stream);
	utility::SerializePrimitive(contig.n_blocks, stream);
	utility::SerializeString(contig.name, stream);
	utility::SerializeString(contig.description, stream);

	size_t size_helper = contig.extra.size();
	utility::SerializePrimitive(size_helper, stream);
	for (uint32_t i = 0; i < contig.extra.size(); ++i) {
		utility::SerializeString(contig.extra[i].first, stream);
		utility::SerializeString(contig.extra[i].second, stream);
	}
	return(stream);
}

std::istream& operator>>(std::istream& stream, YonContig& contig) {
	utility::DeserializePrimitive(contig.idx, stream);
	utility::DeserializePrimitive(contig.n_bases, stream);
	utility::DeserializePrimitive(contig.n_blocks, stream);
	utility::DeserializeString(contig.name, stream);
	utility::DeserializeString(contig.description, stream);

	size_t l_extra;
	utility::DeserializePrimitive(l_extra, stream);
	contig.extra.resize(l_extra);
	for (uint32_t i = 0; i < contig.extra.size(); ++i) {
		utility::DeserializeString(contig.extra[i].first, stream);
		utility::DeserializeString(contig.extra[i].second, stream);
	}
	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonContig& contig) {
	SerializePrimitive(contig.idx, buffer);
	SerializePrimitive(contig.n_bases, buffer);
	SerializePrimitive(contig.n_blocks, buffer);
	SerializeString(contig.name, buffer);
	SerializeString(contig.description, buffer);

	size_t size_helper = contig.extra.size();
	SerializePrimitive(size_helper, buffer);
	for (uint32_t i = 0; i < contig.extra.size(); ++i) {
		SerializeString(contig.extra[i].first, buffer);
		SerializeString(contig.extra[i].second, buffer);
	}
	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, YonContig& contig) {
	DeserializePrimitive(contig.idx, buffer);
	DeserializePrimitive(contig.n_bases, buffer);
	DeserializePrimitive(contig.n_blocks, buffer);
	DeserializeString(contig.name, buffer);
	DeserializeString(contig.description, buffer);

	size_t l_extra;
	DeserializePrimitive(l_extra, buffer);
	contig.extra.resize(l_extra);
	for (uint32_t i = 0; i < contig.extra.size(); ++i) {
		DeserializeString(contig.extra[i].first, buffer);
		DeserializeString(contig.extra[i].second, buffer);
	}
	return(buffer);
}


YonInfo::YonInfo() : yon_type(YON_VCF_HEADER_FLAG) {}
YonInfo::YonInfo(const VcfInfo& vcf_info) : VcfInfo(vcf_info), yon_type(YON_VCF_HEADER_FLAG) {
	this->EvaluateType();
}

bool YonInfo::EvaluateType(void) {
	if (this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
	else if (this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
	else if (this->type == "Flag") this->yon_type = YON_VCF_HEADER_FLAG;
	else if (this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
	else if (this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
	else {
		std::cerr << "Illegal header type: " << this->type << std::endl;
		return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& stream, const YonInfo& info) {
	utility::SerializePrimitive(info.idx, stream);
	utility::SerializeString(info.id, stream);
	utility::SerializeString(info.number, stream);
	utility::SerializeString(info.type, stream);
	utility::SerializeString(info.description, stream);
	utility::SerializeString(info.source, stream);
	utility::SerializeString(info.version, stream);

	return(stream);
}

std::istream& operator>>(std::istream& stream, YonInfo& info) {
	utility::DeserializePrimitive(info.idx, stream);
	utility::DeserializeString(info.id, stream);
	utility::DeserializeString(info.number, stream);
	utility::DeserializeString(info.type, stream);
	utility::DeserializeString(info.description, stream);
	utility::DeserializeString(info.source, stream);
	utility::DeserializeString(info.version, stream);
	info.EvaluateType();

	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonInfo& info) {
	SerializePrimitive(info.idx, buffer);
	SerializeString(info.id, buffer);
	SerializeString(info.number, buffer);
	SerializeString(info.type, buffer);
	SerializeString(info.description, buffer);
	SerializeString(info.source, buffer);
	SerializeString(info.version, buffer);

	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, YonInfo& info) {
	DeserializePrimitive(info.idx, buffer);
	DeserializeString(info.id, buffer);
	DeserializeString(info.number, buffer);
	DeserializeString(info.type, buffer);
	DeserializeString(info.description, buffer);
	DeserializeString(info.source, buffer);
	DeserializeString(info.version, buffer);
	info.EvaluateType();

	return(buffer);
}

YonFormat::YonFormat() : yon_type(YON_VCF_HEADER_FLAG) {}
YonFormat::YonFormat(const VcfFormat& vcf_format) : VcfFormat(vcf_format), yon_type(YON_VCF_HEADER_FLAG) {
	this->EvaluateType();
}

bool YonFormat::EvaluateType(void) {
	if (this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
	else if (this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
	else if (this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
	else if (this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
	else {
		std::cerr << "Illegal header type: " << this->type << std::endl;
		return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& stream, const YonFormat& fmt) {
	utility::SerializePrimitive(fmt.idx, stream);
	utility::SerializeString(fmt.id, stream);
	utility::SerializeString(fmt.number, stream);
	utility::SerializeString(fmt.type, stream);
	utility::SerializeString(fmt.description, stream);

	return(stream);
}

std::istream& operator>>(std::istream& stream, YonFormat& fmt) {
	utility::DeserializePrimitive(fmt.idx, stream);
	utility::DeserializeString(fmt.id, stream);
	utility::DeserializeString(fmt.number, stream);
	utility::DeserializeString(fmt.type, stream);
	utility::DeserializeString(fmt.description, stream);
	fmt.EvaluateType();

	return(stream);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const YonFormat& fmt) {
	SerializePrimitive(fmt.idx, buffer);
	SerializeString(fmt.id, buffer);
	SerializeString(fmt.number, buffer);
	SerializeString(fmt.type, buffer);
	SerializeString(fmt.description, buffer);

	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, YonFormat& fmt) {
	DeserializePrimitive(fmt.idx, buffer);
	DeserializeString(fmt.id, buffer);
	DeserializeString(fmt.number, buffer);
	DeserializeString(fmt.type, buffer);
	DeserializeString(fmt.description, buffer);
	fmt.EvaluateType();

	return(buffer);
}

namespace utility {

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const uint8_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble((uint32_t)data[0]);
	for (uint32_t i = 1; i < n_data; ++i) {
		buffer += ',';
		buffer.AddReadble((uint32_t)data[i]);
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const uint16_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for (uint32_t i = 1; i < n_data; ++i) {
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const uint32_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for (uint32_t i = 1; i < n_data; ++i) {
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const uint64_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	buffer.AddReadble(data[0]);
	for (uint32_t i = 1; i < n_data; ++i) {
		buffer += ',';
		buffer.AddReadble(data[i]);
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const int8_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if (data[0] == YON_BYTE_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (data[0] == YON_BYTE_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (data[i] == YON_BYTE_MISSING) buffer += ",.";
		else if (data[i] == YON_BYTE_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const int16_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}
	// If the first value is end-of-vector then return
	if (data[0] == YON_SHORT_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (data[0] == YON_SHORT_MISSING) buffer += '.';
	else buffer.AddReadble((int32_t)data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (data[i] == YON_SHORT_MISSING) buffer += ",.";
		else if (data[i] == YON_SHORT_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble((int32_t)data[i]);
		}
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const int32_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if (data[0] == YON_INT_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (data[0] == YON_INT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (data[i] == YON_INT_MISSING) buffer += ",.";
		else if (data[i] == YON_INT_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const int64_t* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	// If the first value is end-of-vector then return
	if (data[0] == YON_INT_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (data[0] == YON_INT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (data[i] == YON_INT_MISSING) buffer += ",.";
		else if (data[i] == YON_INT_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

// Special case
yon_buffer_t& ToVcfString_char(yon_buffer_t& buffer, const char* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	buffer += data[0];
	for (uint32_t i = 1; i < n_data; ++i) {
		buffer += ',';
		buffer += data[i];
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const float* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	const uint32_t* const ref = reinterpret_cast<const uint32_t* const>(data);

	// If the first value is end-of-vector then return
	if (ref[0] == YON_FLOAT_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (ref[0] == YON_FLOAT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (ref[i] == YON_FLOAT_MISSING) buffer += ",.";
		else if (ref[i] == YON_FLOAT_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

yon_buffer_t& ToVcfString(yon_buffer_t& buffer, const double* data, const size_t n_data) {
	if (n_data == 0) {
		buffer += '.';
		return(buffer);
	}

	const uint32_t* const ref = reinterpret_cast<const uint32_t* const>(data);

	// If the first value is end-of-vector then return
	if (ref[0] == YON_FLOAT_EOV) {
		buffer += '.';
		return(buffer);
	}

	// First value
	if (ref[0] == YON_FLOAT_MISSING) buffer += '.';
	else buffer.AddReadble(data[0]);

	// Remainder values
	for (uint32_t i = 1; i < n_data; ++i) {
		if (ref[i] == YON_FLOAT_MISSING) buffer += ",.";
		else if (ref[i] == YON_FLOAT_EOV) { return buffer; }
		else {
			buffer += ',';
			buffer.AddReadble(data[i]);
		}
	}

	return(buffer);
}

int32_t* FormatDataHtslib(const uint8_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint16_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint32_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const uint64_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const int8_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) {
		if (src[i] == INT8_MIN)  { dst[i] = bcf_int32_missing; }
		else if (src[i] == INT8_MIN+1) { dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int16_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) {
		if (src[i] == INT16_MIN)  { dst[i] = bcf_int32_missing; }
		else if (src[i] == INT16_MIN+1) { dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int32_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) {
		if (src[i] == INT32_MIN)  { dst[i] = bcf_int32_missing; }
		else if (src[i] == INT32_MIN+1) { dst[i] = bcf_int32_vector_end; }
		else dst[i] = src[i];
	}
	return(dst);
}

int32_t* FormatDataHtslib(const int64_t* const src, int32_t* dst, const size_t n_entries) {
	for (uint32_t i = 0; i < n_entries; ++i) dst[i] = src[i];
	return(dst);
}

int32_t* FormatDataHtslib(const float* const src, int32_t* dst, const size_t n_entries) {
	std::cerr << utility::timestamp("ERROR") << "Cannot convert float into integer." << std::endl;
	return(dst);
}

int32_t* FormatDataHtslib(const double* const src, int32_t* dst, const size_t n_entries) {
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
